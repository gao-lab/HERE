library("data.table")
library("stringr")

input.merged_vcf_gz_filename <- snakemake@input[["merged_vcf_gz_filename"]]
input.merged_variant_only_snpEff_vcf_gz_filename <- snakemake@input[["merged_variant_only_snpEff_vcf_gz_filename"]]


output.merged_variant_only_snpEff_reformatted_dt_filename <- snakemake@output[["merged_variant_only_snpEff_reformatted_dt_filename"]]
output.merged_variant_only_snpEff_ANN_split_dt_filename <- snakemake@output[["merged_variant_only_snpEff_ANN_split_dt_filename"]]
output.merged_variant_only_snpEff_ANN_single_match_dt_filename <- snakemake@output[["merged_variant_only_snpEff_ANN_single_match_dt_filename"]]
output.merged_variant_only_snpEff_event_summary_dt_filename <- snakemake@output[["merged_variant_only_snpEff_event_summary_dt_filename"]]
output.merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_filename <- snakemake@output[["merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_filename"]]
output.merged_vcf_reformatted_single_event_only_melt_dt_filename <- snakemake@output[["merged_vcf_reformatted_single_event_only_melt_dt_filename"]]
output.merged_vcf_reformatted_single_event_only_melt_with_snpEff_ANN_split_annotation_dt_filename <- snakemake@output[["merged_vcf_reformatted_single_event_only_melt_with_snpEff_ANN_split_annotation_dt_filename"]]
output.AEI_per_sample_dt_filename <- snakemake@output[["AEI_per_sample_dt_filename"]]

merged.vcf.reformatted.dt <- setnames(fread(cmd=paste(sep="", "zcat ", input.merged_vcf_gz_filename, " | grep -v '^##'"), header=TRUE, na.strings=".:.:.", drop=c("QUAL", "FILTER", "FORMAT")), "#CHROM", "CHROM")[, SUBSET:=sub(pattern="SUBSET=", replacement="", x=INFO)][, INFO:=NULL]

if (nrow(merged.vcf.reformatted.dt) == 0 ){
    ## no record, all files are not useful
    file.create(output.merged_variant_only_snpEff_reformatted_dt_filename)
    file.create(output.merged_variant_only_snpEff_ANN_split_dt_filename)
    file.create(output.merged_variant_only_snpEff_ANN_single_match_dt_filename)
    file.create(output.merged_variant_only_snpEff_event_summary_dt_filename)
    file.create(output.merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_filename)
    file.create(output.merged_vcf_reformatted_single_event_only_melt_dt_filename)
    file.create(output.merged_vcf_reformatted_single_event_only_melt_with_snpEff_ANN_split_annotation_dt_filename)
    file.create(output.AEI_per_sample_dt_filename)
    q(save="no", status=0)
}

merged.variant.only.snpEff.vcf.reformatted.dt <- setnames(fread(cmd=paste(sep="", "zcat ", input.merged_variant_only_snpEff_vcf_gz_filename, " | grep -v '^##'"), header=TRUE, drop=c("QUAL", "FILTER")), "#CHROM", "CHROM")[, `:=`(SUBSET=sub(pattern=".*SUBSET=([^;]+);.*", replacement="\\1", x=INFO), ANN=sub(pattern=".*ANN=([^;]+)($|;.*)", replacement="\\1", x=INFO))][grepl("LOF=", INFO) == TRUE, `:=`(LOF=sub(pattern=".*LOF=([^;]+)($|;.*)", replacement="\\1", x=INFO))][, INFO:=NULL] ## ~ 7 minutes

fwrite(merged.variant.only.snpEff.vcf.reformatted.dt, output.merged_variant_only_snpEff_reformatted_dt_filename)


merged.variant.only.snpEff.vcf.reformatted.ANN.split.dt <- merged.variant.only.snpEff.vcf.reformatted.dt[, list(ANN.single.match=str_split(string=ANN, pattern=",")[[1]]), list(CHROM, POS, ID, REF, ALT, ANN)] ## ~ 3 minutes


fwrite(merged.variant.only.snpEff.vcf.reformatted.ANN.split.dt, output.merged_variant_only_snpEff_ANN_split_dt_filename)


ANN.single.match.split.matrix <- do.call(rbind, str_split(string=merged.variant.only.snpEff.vcf.reformatted.ANN.split.dt[, ANN.single.match], pattern="\\|"))
merged.variant.only.snpEff.vcf.reformatted.ANN.split.annotation.dt <- setnames(data.table(merged.variant.only.snpEff.vcf.reformatted.ANN.split.dt[, list(CHROM, POS, ID, REF, ALT)], ANN.single.match.split.matrix), 6:21, c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO"))


merged.variant.only.snpEff.vcf.reformatted.ANN.split.annotation.dt[Feature_Type == "transcript", event:=sub(pattern=".*[0-9]+", replacement="\\1", x=HGVS.c)]

## merged.variant.only.snpEff.vcf.reformatted.ANN.split.annotation.dt <- merged.variant.only.snpEff.vcf.reformatted.ANN.split.dt[, data.table(.SD, fread(text=ANN.single.match, sep="|", header=FALSE, col.names=c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO")))]

merged.variant.only.snpEff.vcf.reformatted.ANN.split.annotation.dt[, event:=sub(pattern=".*[0-9]+", replacement="", x=HGVS.c)]

fwrite(merged.variant.only.snpEff.vcf.reformatted.ANN.split.annotation.dt, output.merged_variant_only_snpEff_ANN_single_match_dt_filename)


merged.variant.only.snpEff.vcf.reformatted.event.summary.dt <- merged.variant.only.snpEff.vcf.reformatted.ANN.split.annotation.dt[, list(event.summary=paste(collapse=";", sort(unique(event)))), list(CHROM, POS, ID, REF, ALT)]

fwrite(merged.variant.only.snpEff.vcf.reformatted.event.summary.dt, output.merged_variant_only_snpEff_event_summary_dt_filename)

## now for the sample-dependent tables

## 1. original and molten table with snpEff annotation

merged.vcf.reformatted.with.snpEff.ANN.split.annotation.dt <- merge(x=merged.vcf.reformatted.dt, y=merged.variant.only.snpEff.vcf.reformatted.ANN.split.annotation.dt, by=c("CHROM", "POS", "ID", "REF", "ALT"), all.x=TRUE, all.y=TRUE)


fwrite(merged.vcf.reformatted.with.snpEff.ANN.split.annotation.dt, output.merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_filename)



merged.vcf.reformatted.single.event.only.melt.dt <- melt.data.table(merged.vcf.reformatted.dt[ALT %in% c("A", "C", "G", "T")], id.vars=c("CHROM", "POS", "ID", "REF", "ALT", "SUBSET"), na.rm=TRUE, variable.name="sample", value.name="AC.AN.AF")[, data.table(CHROM, POS, ID, REF, ALT, SUBSET, sample, setnames(data.table(str_split(string=AC.AN.AF, pattern=":", simplify=TRUE)), c("AC", "AN", "AF")))][, `:=`(AC=as.integer(AC), AN=as.integer(AN), AF=as.numeric(AF))]

fwrite(merged.vcf.reformatted.single.event.only.melt.dt, output.merged_vcf_reformatted_single_event_only_melt_dt_filename)




merged.vcf.reformatted.single.event.only.melt.with.snpEff.ANN.split.annotation.dt <- merge(x=merged.vcf.reformatted.single.event.only.melt.dt, y=merged.variant.only.snpEff.vcf.reformatted.ANN.split.annotation.dt, by=c("CHROM", "POS", "ID", "REF", "ALT"), all.x=TRUE, all.y=FALSE, allow.cartesian=TRUE)

fwrite(merged.vcf.reformatted.single.event.only.melt.with.snpEff.ANN.split.annotation.dt, output.merged_vcf_reformatted_single_event_only_melt_with_snpEff_ANN_split_annotation_dt_filename)

## 2. per-sample A-to-G AEI

merged.vcf.reformatted.A.to.G.sites.only.dt <- merged.vcf.reformatted.dt[ID %in% merged.variant.only.snpEff.vcf.reformatted.event.summary.dt[event.summary=='A>G', ID]]


merged.vcf.reformatted.AEI.per.sample.dt <- melt(merged.vcf.reformatted.A.to.G.sites.only.dt[, lapply(.SD, function(temp.column.vector){return(sum(is.na(temp.column.vector) == FALSE)/.N)}), .SDcols=setdiff(colnames(merged.vcf.reformatted.A.to.G.sites.only.dt), c("CHROM", "POS", "ID", "REF", "ALT", "SUBSET"))], measure.vars=setdiff(colnames(merged.vcf.reformatted.A.to.G.sites.only.dt), c("CHROM", "POS", "ID", "REF", "ALT", "SUBSET")), variable.name="sample", value.name="AEI")

fwrite(merged.vcf.reformatted.AEI.per.sample.dt, output.AEI_per_sample_dt_filename)
