library(data.table)
library(glue)

snakemake@input[["stats"]] -> stats.filename
snakemake@input[["coverage"]] -> coverage.filename
snakemake@input[["merged_trim_summary"]] -> merged.trim.summary.filename

if (FALSE) {
    "result/S90_1__check_recal_bam_stat/paired-125-125/200902-GSE101571-full-125-125/GSM2706233/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/120/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/stats" -> stats.filename
    "result/S90_1__check_recal_bam_stat/paired-125-125/200902-GSE101571-full-125-125/GSM2706233/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/120/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/coverage" -> coverage.filename
    "result/S91_1__merge_trim_summary/paired-125-125/200902-GSE101571-full-125-125/GSM2706233/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/merged.trim.summary" -> merged.trim.summary.filename
}

merged.trim.summary.filename -> .;
readLines(.) -> .;
## grep("BEGINNING_OF_SUMMARY_FILE:.*r(1|2|)\\.fastq\\.gz_trimming_report\\.txt", ., value=TRUE) -> trim.summary.filename.lines
grep("Reads written \\(passing filters\\):", ., value=TRUE) -> all.reads.written.lines
all.reads.written.lines -> .;
## all.reads.written.lines[grepl("r(1|)\\.fastq\\.gz_trimming_report\\.txt", trim.summary.filename.lines)] -> all.per.run.read.written.lines
## 
## all.per.run.read.written.lines -> .;
sub(pattern="Reads written \\(passing filters\\):\\s*([0-9,]+)\\s*\\(.*", replacement="\\1", x=.) -> .;
gsub(pattern=",", replacement="", x=.) -> .;
## sum all as integers
as.integer(.) -> .;
sum(.) -> .;
. -> count.of.trimmed.reads


stats.filename -> .;
readLines(.) -> .;
grep("reads mapped:", ., value=TRUE) -> .;
sub(pattern=".*reads mapped:\t([0-9]+)", replacement="\\1", x=.) -> .;
as.integer(.) -> .;
. -> count.of.mapped.reads

coverage.filename -> .;
fread(.) -> .;
setnames(., "#rname", "rname") -> .;
.[rname %in% glue("chr{c(1:22, 'X', 'Y')}"), ] -> .;
.[, region_len := endpos - startpos + 1] -> .;
.[, summed_coverage := meandepth * region_len] -> .;
.[, sum(summed_coverage) / sum(region_len)] -> .;
##
. -> mean.depth.across.whole.genome

data.table(
    count.of.trimmed.reads = count.of.trimmed.reads,
    count.of.mapped.reads = count.of.mapped.reads,
    mean.depth.across.whole.genome = mean.depth.across.whole.genome
) -> .;
.[, mapping.rate:=count.of.mapped.reads / count.of.trimmed.reads] -> .;
## add the rest meta
.[, TYPE:=snakemake@wildcards[["TYPE"]]] -> .;
.[, DATASET:=snakemake@wildcards[["DATASET"]]] -> .;
.[, SAMPLE:=snakemake@wildcards[["SAMPLE"]]] -> .;
##
. -> summary.dt

fwrite(summary.dt, snakemake@output[["summary_dt"]])
