library("Biostrings")
library("GenomicRanges")

mergeDup <- function(bed, ref, fa) {
  df <- read.table(bed, header = F, stringsAsFactors = F)
  names(df) <- header <- c('chr','start','end','id','score','strand')
  gr <- with(df, GRanges(chr, IRanges(start + 1, end), id = id, score = score, strand = strand))
  reference = readDNAStringSet(ref, format = "fasta")[[1]]
  gr.p <- gr[strand(gr) == "+"]
  gr.n <- gr[strand(gr) == "-"]
  dst.p <- DNAStringSet(DNAStringSet(reference, start = start(gr.p), end = end(gr.p)))
  dst.n <- reverseComplement(DNAStringSet(reference, start = start(gr.n), end = end(gr.n)))
  dst <- unique(c(dst.p, dst.n))
  names(dst) <- paste("FC", 1:length(dst), sep = "")
  writeXStringSet(dst, file = fa, format = "fasta")
  return(dst)
}

write2bed <- function(seqname, matches, strand, file = "", append = FALSE) {
  if (file.exists(file) && !append) 
    message("Existing file ", file, " has been overwritten with 'append=FALSE'") 
  hits <- data.frame(seqname = rep.int(seqname, length(matches)), 
                     start = start(matches) - 1, 
                     end = end(matches), 
                     patternID = names(matches), 
                     score = rep.int(0, length(matches)), 
                     strand = rep.int(strand, length(matches)), 
                     check.names = FALSE)
  write.table(hits, 
              file = file, 
              append = append, 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE)
}

queryPattern <- function(query, subject, seqname, tofile = "") {
  # Finds all matches of a pattern in a reference sequence
  cat(">>> Finding all hits in", seqname, "\n")
  ref = readDNAStringSet(subject, format = "fasta")[[1]]
  append <- FALSE
  for (i in seq_len(length(query))) {
    patternID <- names(query)[i]
    pattern <- query[[i]]
    plus.matches <- matchPattern(pattern, ref)
    names(plus.matches) <- rep.int(patternID, length(plus.matches))
    write2bed(seqname = seqname, 
              plus.matches, "+", 
              file = tofile, 
              append = append)
    
    append <- TRUE
    rcpattern <- reverseComplement(pattern)
    minus.matches <- matchPattern(rcpattern, ref)
    names(minus.matches) <- rep.int(patternID, length(minus.matches))
    write2bed(seqname = seqname, 
              minus.matches, "-", 
              file = tofile, 
              append = append)
  }
  cat(">>> DONE! \n")
}

generate_table <- function(bed, fa, out) {
  b <- read.table(bed, header = F, stringsAsFactors = F)
  f <- read.table(fa, header = F, stringsAsFactors = F)
  b$V7 = f$V2[match(b$V4, f$V1)]
  b$V8 = sapply(b$V7, nchar)
  write.table(b, file = out, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}


#############################################################################################################
#setwd("~/cosRNA/")
#cp <- mergeDup(bed = "cp.bed", ref = "ref/TAIR10_chrC.fas", fa = "cp_uniq.fa")
#queryPattern(query = cp, subject = "ref/TAIR10_chrC.fas", seqname = "chloroplast", tofile = "cp_uniq.bed")

#mt <- mergeDup(bed = "mt.bed", ref = "ref/TAIR10_chrM.fas", fa = "mt_uniq.fa")
#queryPattern(query = mt, subject = "ref/TAIR10_chrM.fas", seqname = "mitochondria", tofile = "mt_uniq.bed")

## convert cp_uniq.fa and mt_uniq.fa to csv format
#generate_table(bed = "cp_uniq.bed", fa = "cp_uniq.txt", out = "chloroplast.txt")
#generate_table(bed = "mt_uniq.bed", fa = "mt_uniq.txt", out = "mitochondria.txt")
