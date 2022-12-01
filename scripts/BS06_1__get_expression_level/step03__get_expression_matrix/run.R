library("ballgown")
library(data.table)
params.ballgown_directory_filenames_collection <- snakemake@params[["ballgown_directory_filenames_collection"]]
params.temp_linking_directory <- snakemake@params[["temp_linking_directory"]]
output.combined_texpr_FPKM_matrix_filelname <- snakemake@output[["combined_texpr_FPKM_matrix_filename"]]
output.combined_gexpr_FPKM_matrix_filelname <- snakemake@output[["combined_gexpr_FPKM_matrix_filename"]]

## print(1)
## print(snakemake@params)
print(params.ballgown_directory_filenames_collection)
print(params.temp_linking_directory)
system(paste(sep="", "rm -fr ", params.temp_linking_directory, "/*"))
dir.create(params.temp_linking_directory, recursive=TRUE)
sample.accessions.vector <- sub(pattern="^result/S06_1__get_expression_level/[^/]+/[^/]+/([^/]+)/__merged__/.*$", replacement="\\1", x=params.ballgown_directory_filenames_collection)




## must name the samples' LAST DIRECTORY correctly before loading them into ballgown
## otherwise the ballgown gexpr and texpr will use the last directory name of each sample path, and err if any two of these paths coincide (in which case only the first path is used twice)

not.used.variable <- mapply(params.ballgown_directory_filenames_collection, sample.accessions.vector, FUN=function(temp.ballgown.directory.filename, temp.sample.accession){
    file.symlink(normalizePath(temp.ballgown.directory.filename), paste(sep="", params.temp_linking_directory, "/", temp.sample.accession))
    return(NULL)
})

ballgown.result <- ballgown(samples=paste(sep="", params.temp_linking_directory, "/", sample.accessions.vector), meas="all")



## get transcript expression matrix
ballgown.result.texpr.FPKM.matrix <- texpr(ballgown.result, meas="FPKM")
rownames(ballgown.result.texpr.FPKM.matrix) <- transcriptNames(ballgown.result)
colnames(ballgown.result.texpr.FPKM.matrix) <- sample.accessions.vector

write.table(ballgown.result.texpr.FPKM.matrix, output.combined_texpr_FPKM_matrix_filelname, row.names=TRUE, col.names=TRUE)

## get gene expression matrix
ballgown.result.gexpr.FPKM.matrix <- gexpr(ballgown.result)
colnames(ballgown.result.gexpr.FPKM.matrix) <- sample.accessions.vector
write.table(ballgown.result.gexpr.FPKM.matrix, output.combined_gexpr_FPKM_matrix_filelname, row.names=TRUE, col.names=TRUE)
