library("Biostrings")
library("stringr")
library("magrittr")
library("foreach")
library("data.table")

contigs.DNAStringSet <- readDNAStringSet(snakemake@input[["contigs_fasta_filename"]])
if (FALSE) {
    contigs.DNAStringSet <- readDNAStringSet("external/contigs/hg38.fa")
}

temp.result.directory <- snakemake@params[["temp_result_directory"]]
if (FALSE){
    temp.result.directory <- "result/sA02_4__generate_all_transcribable_Alu_A_vcf/hg38.fa/temp/"
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

## annotate with snpEff
snpEff.config.filename <- snakemake@params[["snpEff_config_filename"]]
if (FALSE){
    snpEff.config.filename <- "result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/hg38.fa/32/snpEff.config"
}

snpEff.database.name <- snakemake@params[["snpEff_database_name"]]
if (FALSE){
   snpEff.database.name <- "hg38.fa.GENCODE.32"
}

snpEff.Xmx <- snakemake@params[["snpEff_Xmx"]]
if (FALSE){
    snpEff.Xmx <- "-Xmx60G"
}

all.contigs.Alu.only.combined.snpEff.annotated.vcf.gz.filename <- snakemake@output[["all_contigs_Alu_only_combined_snpEff_annotated_vcf_gz_filename"]]
if (FALSE) {
    all.contigs.Alu.only.combined.snpEff.annotated.vcf.gz.filename <- "result/sA02_4__generate_all_transcribable_Alu_A_vcf/hg38.fa/all.contigs.Alu.only.combined.snpEff.annotated.vcf.gz"
}

system(paste(sep="", "snpEff ann ", snpEff.Xmx, " -verbose -config ", snpEff.config.filename, " ", snpEff.database.name, " ", all.contigs.Alu.only.combined.vcf.gz.filename, " | bcftools view --threads ", threads, " --output-file ", all.contigs.Alu.only.combined.snpEff.annotated.vcf.gz.filename, " -Oz "))

## get reformatted table
all.contigs.Alu.only.combined.snpEff.annotated.vcf.reformatted.dt <- fread(cmd=paste(sep="", "zcat ", all.contigs.Alu.only.combined.snpEff.annotated.vcf.gz.filename, " | grep -v '^##'"), header=TRUE, drop=c("QUAL", "FILTER")) %>%
    setnames("#CHROM", "CHROM") %>%
    {.[, `:=`(ANN=sub(pattern=".*ANN=([^;]+)($|;.*)", replacement="\\1", x=INFO))]} %>%
    {.[, list(ANN.single.match=str_split(string=ANN, pattern=",")[[1]]), list(CHROM, POS, ID, REF, ALT, ANN)]} %>%
    {
        ANN.single.match.split.matrix <- do.call(rbind, str_split(string=.[, ANN.single.match], pattern="\\|"))
        data.table(.[, list(CHROM, POS, ID, REF, ALT)], ANN.single.match.split.matrix) %>%
            setnames(6:21, c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO"))
    }
