library("data.table")



## first row, sample profile (sequenced only; could have no edits)

phenotype.output.at.gsm.level.dt <- fread(snakemake@input[["phenotype_output_at_gsm_level_dt_filename"]])
if (FALSE) { 
    phenotype.output.at.gsm.level.dt <- fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt")
}

sample.sequenced.dt <- rbindlist(lapply(readLines(paste(sep="", "external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/", snakemake@wildcards[['DATASET_COLLECTION_NAME']])), function(temp.DATASET_RNA_EDITING_NAME){return(unique(fread(paste(sep="", "external/DATASET_RNA_EDITING_NAME_DIRECTORY/", temp.DATASET_RNA_EDITING_NAME))))}))

if (FALSE) {
    sample.sequenced.dt <- rbindlist(lapply(readLines(paste(sep="", "external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/", "210215-sixth-dataset")), function(temp.DATASET_RNA_EDITING_NAME){return(unique(fread(paste(sep="", "external/DATASET_RNA_EDITING_NAME_DIRECTORY/", temp.DATASET_RNA_EDITING_NAME))))}))
}

phenotype.output.at.gsm.level.sequenced.only.dt <- phenotype.output.at.gsm.level.dt[gsm %in% sample.sequenced.dt[, SAMPLE_NAME]]


total.sample.count.for.normal.stages.dt <- phenotype.output.at.gsm.level.sequenced.only.dt[is.normal==TRUE, list(total.sample.count=.N), list(stage)]

fwrite(total.sample.count.for.normal.stages.dt, snakemake@output[["total_sample_count_for_normal_stages_dt_csv_filename"]])
if (FALSE) {
    fwrite(total.sample.count.for.normal.stages.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.sample.count.for.normal.stages.dt.csv")
}

total.sample.count.for.GSE133854.all.dt <- phenotype.output.at.gsm.level.sequenced.only.dt[gse == "GSE133854", list(total.sample.count=.N), list(stage, is.normal, disease)][is.normal==TRUE, disease:="biparental"]

fwrite(total.sample.count.for.GSE133854.all.dt, snakemake@output[["total_sample_count_for_GSE133854_all_dt_csv_filename"]])
if (FALSE) {
    fwrite(total.sample.count.for.GSE133854.all.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.sample.count.for.GSE133854.all.dt.csv")
}




