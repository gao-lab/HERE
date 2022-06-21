library("data.table")
library("stringr")
library("magrittr")

pileup.filename <- snakemake@input[["alignment_pileup_filename"]]

if (FALSE){
    pileup.filename <- "result/PS15_1__get_sample_RNA_editing_sites_v3/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/GSM1160112/__merged__/filter-by-soapnuke/hg19.fa/Ensembl75.GRCh37/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-3.6.0/GATK-3.6.0-default/human.b144.GRCh37p13/common_all/alignment.pileup.gz"
}

pileup.dt <- fread(pileup.filename, header=FALSE, col.names=c("chr", "pos", "ref", "count", "bases", "qualities"))
date() ## 16:01-16:06

## remove the following to keep _mapped_ bases only: 
## 1. "\\^." (start of read, and its mapping quality)
## 2. "\\+([0-9]+).{n}" (insertion from the read after this base with length n == \1 )
## 3. "-([0-9]+).{n}" (deletion on the read after this base with length n == \1 )
## 4. "\\*" (deletion on the current base; no need to filter against "#" because we did not use "--reverse-del" in samtools mpileup)
## 5. "[><]" (cigar N)
## 6. "$" (last position of a read)
## 


date(); pileup.dt[, contains.indel:=grepl("[-+]", bases)]; date()

## test set for "-"
## temp.str <- pileup.dt[c(112988,303502,875725,878187), bases]

{
    str_locate_all(string=temp.str, pattern="-[0-9]+")
    str_extract_all(string=temp.str, pattern="-[0-9]+") %>% {lapply(., FUN=sub, pattern="^-", replacement="")}
}

pileup.dt[, ]
