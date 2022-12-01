library("data.table")
library("magrittr")

input.gencode.3utr.dt.csv.gz.filename <- snakemake@input[["gencode_3utr_dt_csv_gz_filename"]]
input.snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt.txt.gz.filename <- snakemake@input[["snpEff_annotation_for_subset_recurrent_edits_on_valid_transcripts_dt_txt_gz_filename"]]
output.gencode.3utr.and.edit.CJ.dt.csv.gz.filename <- snakemake@output[["gencode_3utr_and_edit_CJ_dt_csv_gz_filename"]]

gencode.3utr.valid.dt <- fread(input.gencode.3utr.dt.csv.gz.filename)

snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt <- fread(input.snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt.txt.gz.filename)

gencode.3utr.and.edit.CJ.dt <- merge(x=gencode.3utr.valid.dt[, list(CHROM, UTR.block.start=start, UTR.block.end=end, UTR.strand=strand, transcript_id)],
      y=snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt[Annotation=='3_prime_UTR_variant'][, list(CHROM, edit.POS=POS, Feature_ID, HGVS.c)],
      by.x=c("CHROM", "transcript_id"),
      by.y=c("CHROM", "Feature_ID"),
      all.x=TRUE, all.y=TRUE, allow.cartesian=TRUE) %>%
    ## mark within-range edits
    {.[, edit.is.in.UTR.block:=(edit.POS >= UTR.block.start & edit.POS <= UTR.block.end)]} %>%
    ## get relative position of edit w.r.t. 3'-UTR
    ## note: 3'-UTR contains the stop codon. The number in HGVS.c is the position of the edit relative to **the first base in the modified 3'-UTR where the first three bases for stop codon are discarded**. Therefore the relative position of edit w.r.t. 3'UTR should be "HGVS.c number + 3"
    {.[edit.is.in.UTR.block == TRUE, HGVS.c.integer := sub(pattern="c\\.\\*([0-9]+)[^0-9]+", replacement="\\1", x=HGVS.c) %>% as.integer]} %>%
    {.[edit.is.in.UTR.block == TRUE, edit.rel.POS.wrt.3UTR := HGVS.c.integer + 3]}


fwrite(gencode.3utr.and.edit.CJ.dt, output.gencode.3utr.and.edit.CJ.dt.csv.gz.filename)
