library(data.table)
library(glue)

{

    snakemake@input[["Editing_file_collection"]] -> .;
    rbindlist(lapply(., FUN=function(temp.Editing.file){
        fread(temp.Editing.file) -> TEMP.dt
        TEMP.dt[, list(A2CEditingIndex, A2GEditingIndex, A2TEditingIndex, C2AEditingIndex, C2GEditingIndex, C2TEditingIndex)] -> TEMP.dt
        TEMP.dt[, filename:=temp.Editing.file] -> TEMP.dt
        TEMP.dt[, SAMPLE:=sub(pattern="result/S80_1__compute_AEI/.*/.*/([^/]+)/__merged__.*", replacement="\\1", x=filename)] -> TEMP.dt
        TEMP.dt
    }), use.names=TRUE) -> .
    fwrite(., snakemake@output[["combined_Editing_file"]])


}
