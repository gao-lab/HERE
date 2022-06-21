library("Biostrings")
library("stringr")
library("magrittr")
library("foreach")
library("data.table")

contigs.DNAStringSet <- readDNAStringSet(snakemake@input[["contigs_fasta_filename"]])
if (FALSE) {
    contigs.DNAStringSet <- readDNAStringSet("external/contigs/hg38.fa")
}

temp.result.directory <- snakemake@params[["result_directory"]]
if (FALSE){
    temp.result.directory <- "result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/"
}

repeatmasker.repFamily.Alu.bed.filename <- snakemake@input[["repeatmasker_repFamily_Alu_bed_filename"]]
if (FALSE){
    repeatmasker.repFamily.Alu.bed.filename <- "external/UCSC.Table.Browser.repeatmasker/repFamily.Alu/repeatmasker.bed"
}

threads <- as.integer(snakemake@threads)
if (FALSE){
    threads <- 20
}

temp.result.directory %>%
    {system(paste(sep="", "rm -fr ", .)); .} %>%
    {dir.create(path=., recursive=TRUE)}


## create vcf
not.used.variable <- foreach(temp.contig.DNAString=contigs.DNAStringSet, temp.contig.name=names(contigs.DNAStringSet)) %do% {
    cat(date(), " Processing contig ", temp.contig.name, "\n")
    temp.contig.A.position <- str_locate_all(string=temp.contig.DNAString, pattern="A")[[1]][, 1]
    temp.contig.T.position <- str_locate_all(string=temp.contig.DNAString, pattern="T")[[1]][, 1]
    temp.vcf.filename <- paste(sep="", temp.result.directory, "/", temp.contig.name, ".vcf")
    writeLines(c(
        "##fileformat=VCFv4.2",
        paste(sep="", "##contig=<ID=", temp.contig.name, ">")
    ), temp.vcf.filename)
    list(
        data.table(pos=temp.contig.A.position, ref="A", alt="G"),
        data.table(pos=temp.contig.T.position, ref="T", alt="C")
    ) %>%
        rbindlist %>%
        {.[order(pos)]} %>%
        {.[, list("#CHROM"=temp.contig.name, POS=pos, ID=".", REF=ref, ALT=alt, QUAL=".", FILTER=".", INFO=".")]} %>%
        {fwrite(x=., file=temp.vcf.filename, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)}
    system(paste(sep="", "bgzip --threads ", threads, " ", temp.vcf.filename))
    NULL
}

## index vcf
not.used.variable <- foreach(temp.contig.name=names(contigs.DNAStringSet)) %do% {
    cat(date(), " Indexing vcf.gz for contig ", temp.contig.name, "\n")
    temp.vcf.gz.filename <- paste(sep="", temp.result.directory, "/", temp.contig.name, ".vcf.gz")
    system(paste(sep="", "tabix ", temp.vcf.gz.filename))
    NULL
}

Sys.setenv("PATH"=paste(sep="", "/gpfs2/gaog_pkuhpc/users/dingy/Project.collaborated/Project.embryo.I/tools/bcftools-1.10.2:", Sys.getenv("PATH")))

## subset Alu and index
not.used.variable <- foreach(temp.contig.name=names(contigs.DNAStringSet)) %do% {
    cat(date(), " Subsetting Alu and indexing for contig ", temp.contig.name, "\n")
    temp.vcf.gz.filename <- paste(sep="", temp.result.directory, "/", temp.contig.name, ".vcf.gz")
    temp.Alu.only.vcf.gz.filename <- paste(sep="", temp.result.directory, "/", temp.contig.name, ".Alu.only.vcf.gz")
    system(paste(sep="", "bcftools view --threads ", threads, " --regions-file ", repeatmasker.repFamily.Alu.bed.filename, " --output-file ", temp.Alu.only.vcf.gz.filename , " -Oz ", temp.vcf.gz.filename))
    system(paste(sep="", "tabix ", temp.Alu.only.vcf.gz.filename))
}

## combine into a single vcf
temp.all.contigs.Alu.only.vcf.gz.filename.list.txt.filename <- paste(sep="", temp.result.directory, "/temp.all.contigs.Alu.only.vcf.gz.filename.list.txt")
all.contigs.Alu.only.combined.vcf.gz.filename <- paste(sep="", temp.result.directory, "/all.contigs.Alu.only.combined.vcf.gz")
names(contigs.DNAStringSet) %>%
    {paste(sep="", temp.result.directory, "/", ., ".Alu.only.vcf.gz")} %>%
    {writeLines(., temp.all.contigs.Alu.only.vcf.gz.filename.list.txt.filename)}
system(paste(sep="", "bcftools concat --threads ", threads, " --file-list ", temp.all.contigs.Alu.only.vcf.gz.filename.list.txt.filename, " --output ", all.contigs.Alu.only.combined.vcf.gz.filename, " -Oz "))
system(paste(sep="", "tabix ", all.contigs.Alu.only.combined.vcf.gz.filename))
