sRNA miner
=========================================

1. In order to use the shiny app,  you have to put bam and index files ( or their soft link) in the **dat** folder. 
> BAM file is required to be sorted by position and indexed. 
>The BAM index file should be named by appending '.bai' to the BAM file name and it must reside in the same directory as the BAM file.

2. It may takes several minutes to get the coverage infomation, **DO NOT REFRESH** during this time.
> If you have a huge bam file, please split by chromosomes or make sure you have enough RAM before running sRNA miner.

3. The resulting annotation file can be downloaded in the **Report** tab in *bed* format.  
> If a sorted bed file is needed, one option is to use external utility `sort`.  

````bash
# sort by chromosome then start position
sort -k 1,1 -k 2n,2 -o sorted.bed unsorted.bed
````
