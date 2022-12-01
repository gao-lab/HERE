library("data.table")
library("magrittr")

input.reference.GTF.filename <- snakemake@input[["reference_GTF_filename"]]
if (FALSE) {
    input.reference.GTF.filename <- "./external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf"
}
output.gencode.3utr.dt.csv.gz.filename <- snakemake@output[["gencode_3utr_dt_csv_gz_filename"]]
if (FALSE) {
    output.gencode.3utr.dt.csv.gz.filename <- "result/A02_8__get_editing_effect_on_miRNA_binding_sites/step01__get_UTR_blocks/32/gencode.3utr.dt.csv.gz"
}
output.gencode.3utr.gff.filename <- snakemake@output[["gencode_3utr_gff_filename"]]
if (FALSE) {
    output.gencode.3utr.gff.filename <- "result/A02_8__get_editing_effect_on_miRNA_binding_sites/step01__get_UTR_blocks/32/gencode.3utr.gff"
}

if (FALSE) {
    print(input.reference.GTF.filename)
    print(output.gencode.3utr.dt.csv.gz.filename)
    print(output.gencode.3utr.gff.filename)
}

input.dt <- fread(cmd=paste(sep="", "grep -v '#' ", input.reference.GTF.filename), select=c(1,3,4,5,7,9), col.names=c("CHROM", "type", "start", "end", "strand", "group"))

## pick core cols
gencode.utr.dt <-  input.dt %>%
    ## get stop_codon and UTR rows
    {.[type %in% c("stop_codon", "UTR")]} %>%
    ## get gene_id, transcript id and remove group
    {.[, gene_id:=sub(pattern=".*gene_id \"([^\"]+)\".*", replacement="\\1", x=group)]} %>%
    {.[, transcript_id:=sub(pattern=".*transcript_id \"([^\"]+)\".*", replacement="\\1", x=group)]} %>%
    {.[, group:=NULL]} %>%
    ## keep only those transcripts with stop codon (see ENST00000602680.1 for an example of stop codon-free transcript)
    {.[transcript_id %in% .[type=="stop_codon", transcript_id] ]} %>%
    ## mark the first 3 prime UTR block for each transcript
    ## NOTE: some stop codons are split by splicing sites, and therefore we need to fetch the leftmost start/rightmost end for each transcript on plus/minus strand
    {.[strand == "+", is.first.3UTR.region:= (type == "UTR") & (start == .SD[type=="stop_codon", min(start)]), transcript_id]} %>%
    {.[strand == "-", is.first.3UTR.region:= (type == "UTR") & (end == .SD[type=="stop_codon", max(end)]), transcript_id]} %>%
    ## mark the whole 3 prime UTR region for each transcript
    {.[type == "UTR" & strand == "+", is.3UTR.region:=(start >= .SD[is.first.3UTR.region==TRUE, start]), transcript_id]} %>%
    {.[type == "UTR" & strand == "-", is.3UTR.region:=(end <= .SD[is.first.3UTR.region==TRUE, end]), transcript_id]} %>%
    ## return results
    identity


gencode.3utr.dt <- gencode.utr.dt[is.3UTR.region == TRUE]

fwrite(gencode.3utr.dt, output.gencode.3utr.dt.csv.gz.filename)

gencode.3utr.gff.dt <- gencode.3utr.dt[, list(
    seqname=CHROM,
    source="GENCODE.derived",
    feature="exon",
    start=start,
    end=end,
    score=".",
    strand=strand,
    frame=".",
    group=paste(sep="", "gene_id \"", gene_id, "\"; transcript_id \"", transcript_id, "\"")
)]


## unit test passed
## fwrite(gencode.3utr.gff.dt[grepl("ENST00000641515.2|ENST00000618323.4|ENST00000417435.5|ENST00000426406.3", group)], "./temp.miRNA/step2.gencode.3utr.gtf", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

fwrite(gencode.3utr.gff.dt, output.gencode.3utr.gff.filename, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
