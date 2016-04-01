# This should detect and install missing packages before loading them.

pkg <- c("shiny", "Rsamtools", "GenomicRanges", "GenomicAlignments", "rmarkdown")
npkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(npkg)) install.packages(npkg)
lapply(pkg, function(x) {
  library(x, character.only = T)
})

getSeqnamesBam <- function(file) {
  bf <- BamFile(file)
  seqnames(seqinfo(bf))
}

getCoverage <- function(bam, chr) {
  if (!file.exists(paste(bam, "bai", sep = "."))) {
    stop("Unable to find index")
  }
  gal <- readGAlignments(bam)
  gr <- sort(granges(gal[seqnames(gal) %in% chr]))
  rm(gal)
  
  len <- seqlengths(gr)[[chr]]
  
  cvg.p5 <- coverage(resize(gr[strand(gr) == "+"], 1, "start"))
  cvg.p3 <- coverage(resize(gr[strand(gr) == "+"], 1, "end"))
  
  cvg.n5 <- coverage(resize(gr[strand(gr) == "-"], 1, "start"))
  cvg.n3 <- coverage(resize(gr[strand(gr) == "-"], 1, "end"))
  
  cvg.p5 <- as.numeric(cvg.p5[[chr]])
  cvg.p3 <- as.numeric(cvg.p3[[chr]])
  
  cvg.n5 <- as.numeric(cvg.n5[[chr]])
  cvg.n3 <- as.numeric(cvg.n3[[chr]])
  
  pos <- data.frame(seqnames = chr, 
                    start = 1:len, 
                    end = 1:len, 
                    strand = "+", 
                    cvg.p5 = cvg.p5, 
                    cvg.p3 = cvg.p3)
  neg <- data.frame(seqnames = chr, 
                    start = 1:len, 
                    end = 1:len, 
                    strand = "-", 
                    cvg.n5 = cvg.n5, 
                    cvg.n3 = cvg.n3)
  
  return(list(pos = pos, neg = neg))
}

localMax <- function(x, bg, w = 15) {
  len <- length(x)
  idx <- which(x >= bg)
  
  n <- 0
  lm <- 0
  
  for (i in idx) {
    if (i > w & i < len - w) {
      max <- which.max(x[(i - w):(i + w)]) + (i - w - 1)
    } else if (i <= w) {
      max <- which.max(x[1:(i + w)])
    } else {
      max <- which.max(x[(i - w):len]) + (i - w - 1)
    }
    
    if (max == i) {
      n <- n + 1
      lm[n] <- max
    }
  }
  return(lm)
}

winSum <- function(idx, cvg, w = 15) {
  len <- length(cvg)
  j <- 0
  sum <- 0
  
  for (i in idx) {
    j <- j + 1
    if (i - w < 0) {
      s <- 1
      e <- i + w
    } else if (i + w > len) {
      s <- i + w
      e <- len
    } else {
      s <- i - w
      e <- i + w
    }
    
    sum[j] <- sum(cvg[s:e])
  }
  return(sum)
}


seft <- function(idx, cvg, w, min.se) {
  s1 <- winSum(idx, cvg, w = 1)
  s2 <- winSum(idx, cvg, w = w)
  return(idx[which(s1/s2 >= min.se)])
}

matchEnds <- function(lm, cvg, forward = TRUE) {
  if (forward) {
    mate.win <- sapply(lm, function(x) (x + 14):(x + 49))
  } else {
    mate.win <- sapply(lm, function(x) (x - 49):(x - 14))
  }
  mw <- 0
  mate <- 0
  for (i in 1:ncol(mate.win)) {
    mw <- mate.win[, i][mate.win[, i] > 0 & mate.win[, i] <= length(cvg)]
    mate[i] <- mw[1] + which.max(cvg[mw]) - 1
  }
  return(mate)
}

getPeaks <- function(dat, chr, w = 15, bg = 40, min.se = 0.5, min.ft = 0.8, Edges = NULL) {
  cvg.p5 <- dat[["pos"]]$cvg.p5
  cvg.p3 <- dat[["pos"]]$cvg.p3
  cvg.n5 <- dat[["neg"]]$cvg.n5
  cvg.n3 <- dat[["neg"]]$cvg.n3
  
  lm.p5 <- localMax(x = cvg.p5, bg = bg, w = w)
  lm.p3 <- localMax(x = cvg.p3, bg = bg, w = w)
  lm.n5 <- localMax(x = cvg.n5, bg = bg, w = w)
  lm.n3 <- localMax(x = cvg.n3, bg = bg, w = w)
  
  p5 <- seft(lm.p5, cvg.p5, w = w, min.se = min.se)
  p3 <- seft(lm.p3, cvg.p3, w = w, min.se = min.se)
  n5 <- seft(lm.n5, cvg.n5, w = w, min.se = min.se)
  n3 <- seft(lm.n3, cvg.n3, w = w, min.se = min.se)

  p5m <- matchEnds(p5, cvg.p3, forward = TRUE)
  p3m <- matchEnds(p3, cvg.p5, forward = FALSE)
  n5m <- matchEnds(n5, cvg.n3, forward = FALSE)
  n3m <- matchEnds(n3, cvg.n5, forward = TRUE)
  
  gr.p5 <- GRanges(seqnames = chr, IRanges(start = p5, end = p5m), strand = "+")
  gr.p3 <- GRanges(seqnames = chr, IRanges(start = p3m, end = p3), strand = "+")
  gr.n5 <- GRanges(seqnames = chr, IRanges(start = n5m, end = n5), strand = "-")
  gr.n3 <- GRanges(seqnames = chr, IRanges(start = n3, end = n3m), strand = "-")
  
  if (Edges == 1) {
    gr <- c(gr.p5, gr.n5)
  } else if (Edges == 2) {
    gr <- c(gr.p3, gr.n3)
  } else if (Edges == 3) {
    gr <- unique(c(gr.p5, gr.p3, gr.n5, gr.n3))
  } else if (Edges == 4) {
    ## merge overlap footprints, in case there are.
    ov.p <- findOverlaps(gr.p5, gr.p3, minoverlap = 15L)
    ov.n <- findOverlaps(gr.n5, gr.n3, minoverlap = 15L)
    gr <- unique(c(gr.p5[queryHits(ov.p)], gr.n5[queryHits(ov.n)]))
  } else {
    stop("Error: Parameter 'ends = ", Edges, "' is not recognized.")
  }
}

# getCoveragesBam(file, selseqlen){
#   seqname <- names(selseqlen)
#   seqlen <- unname(selseqlen)
#   param <- ScanBamParam(what = c("pos", "qwidth", "strand"),
#                         which = GRanges(seqname, IRanges(1,seqlen)),
#                         flag = scanBamFlag(isUnmappedQuery = FALSE))
#   x <- scanBam(file,param = param)[[1]]
#   coverage(IRanges(x[["pos"]], width = x[["qwidth"]]))
# }
