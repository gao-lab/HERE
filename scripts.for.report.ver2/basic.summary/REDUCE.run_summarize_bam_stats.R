library("data.table")

snakemake@input[["summary_dts_collection"]] -> .;
lapply(., FUN=fread) -> .;
rbindlist(., use.names=TRUE) -> .;
fwrite(., snakemake@output[["combined_summary_dt"]])


