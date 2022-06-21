library("ballgown")
library(data.table)
params.ballgown_directory_filenames_collection <- snakemake@params[["ballgown_directory_filenames_collection"]]
params.temp_linking_directory <- snakemake@params[["temp_linking_directory"]]
output.combined_texpr_cov_matrix_filelname <- snakemake@output[["combined_texpr_cov_matrix_filename"]]
output.combined_gexpr_cov_matrix_filelname <- snakemake@output[["combined_gexpr_cov_matrix_filename"]]

if (FALSE){
    params.ballgown_directory_filenames_collection <- c(
        "result/S06_1__get_expression_level/paired-101-101/201031-GSE130289-full139-101-101/GSM3735310/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/100/STAR-expression/default/stringtie/none/",
        "result/S06_1__get_expression_level/paired-101-101/201031-GSE130289-full139-101-101/GSM3735313/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/100/STAR-expression/default/stringtie/none/",
        "result/S06_1__get_expression_level/paired-101-101/201031-GSE130289-full139-101-101/GSM3735316/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/100/STAR-expression/default/stringtie/none/"
    )
    params.temp_linking_directory <- "result/BS06_1__get_expression_level/210215-sixth-datasets/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/step04__temp_linking/"
}

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
ballgown.result.texpr.cov.matrix <- texpr(ballgown.result, meas="cov")
rownames(ballgown.result.texpr.cov.matrix) <- transcriptNames(ballgown.result)
colnames(ballgown.result.texpr.cov.matrix) <- sample.accessions.vector

write.table(ballgown.result.texpr.cov.matrix, output.combined_texpr_cov_matrix_filelname, row.names=TRUE, col.names=TRUE)

rm(ballgown.result.texpr.cov.matrix)
gc()

## get gene expression matrix
ballgown.result.gexpr.cov.matrix <- gexpr(ballgown.result)
colnames(ballgown.result.gexpr.cov.matrix) <- sample.accessions.vector
write.table(ballgown.result.gexpr.cov.matrix, output.combined_gexpr_cov_matrix_filelname, row.names=TRUE, col.names=TRUE)
