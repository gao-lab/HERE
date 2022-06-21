library("fastqcr")
library("magrittr")
library("foreach")
library("data.table")

fastqc.result.directories.vector <- snakemake@params[["fastqc_result_directory_collection"]]



if (FALSE){
    fastqc.result.directories.vector <- c("result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706237/SRR5837350/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706237/SRR5837349/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706266/SRR5837401/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706266/SRR5837400/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706265/SRR5837399/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706265/SRR5837398/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706245/SRR5837362/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706245/SRR5837361/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706261/SRR5837389/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706261/SRR5837388/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706241/SRR5837356/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706241/SRR5837355/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706244/SRR5837359/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706244/SRR5837360/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706260/SRR5837386/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706260/SRR5837387/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706238/SRR5837351/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706238/SRR5837352/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706262/SRR5837391/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706262/SRR5837390/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706264/SRR5837397/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706264/SRR5837396/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706263/SRR5837393/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706263/SRR5837392/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706263/SRR5837395/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706263/SRR5837394/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706267/SRR5837403/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-100-100/200902-GSE101571-full-100-100/GSM2706267/SRR5837402/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-150-150/200902-GSE101571-full-150-150/GSM3081011/SRR6940545/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-150-150/200902-GSE101571-full-150-150/GSM3081010/SRR6940544/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-125-125/200902-GSE101571-full-125-125/GSM2706235/SRR5837347/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-125-125/200902-GSE101571-full-125-125/GSM2706234/SRR5837346/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-125-125/200902-GSE101571-full-125-125/GSM2706242/SRR5837357/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-125-125/200902-GSE101571-full-125-125/GSM2706243/SRR5837358/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-125-125/200902-GSE101571-full-125-125/GSM2706240/SRR5837354/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-125-125/200902-GSE101571-full-125-125/GSM2706236/SRR5837348/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-125-125/200902-GSE101571-full-125-125/GSM2706239/SRR5837353/default/fastqc_result",
"result/S01_1__check_RNA_Seq/paired-125-125/200902-GSE101571-full-125-125/GSM2706233/SRR5837345/default/fastqc_result")
}

fastqc.result.zip.filenames.vector <- foreach(temp.result.directory=fastqc.result.directories.vector, .combine=c) %do% {
    if (grepl("result/S01_1__check_RNA_Seq/paired", temp.result.directory) == TRUE){
        paste(sep="", temp.result.directory, c("/r1/r1_fastqc.zip", "/r2/r2_fastqc.zip"))
    } else if (grepl("result/S01_1__check_RNA_Seq/single", temp.result.directory) == TRUE){
        paste(sep="", temp.result.directory, "/r/r_fastqc.zip")
    } else {
        stop(paste0("Unsupported type of the directory ", temp.result.directory))
    }
}

fastqc.results.adapter.content.dt <- foreach(temp.zip.filename=fastqc.result.zip.filenames.vector) %do% {
    qc_read(temp.zip.filename, modules='Adapter Content')[['adapter_content']] %>% data.table %>% {.[, zip.filename:=temp.zip.filename]}
} %>% rbindlist

fwrite(fastqc.results.adapter.content.dt, snakemake@output[["adapter_content_dt_txt_filename"]])


fastqc.results.overrepresented.sequences.dt <- foreach(temp.zip.filename=fastqc.result.zip.filenames.vector) %do% {
    qc_read(temp.zip.filename, module='Overrepresented sequences')[['overrepresented_sequences']] %>% {if (ncol(.)==4) {data.table(.)} else {data.table()}} %>% {.[, zip.filename:=temp.zip.filename]}
} %>% rbindlist(fill=TRUE, use.names=TRUE)

fwrite(fastqc.results.overrepresented.sequences.dt, snakemake@output[["overrepresented_sequences_dt_txt_filename"]])
