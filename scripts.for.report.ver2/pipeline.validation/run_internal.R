library("data.table")
library("ggpubr")
library("foreach")
library("iterators")
source("./scripts/common/ggpubr.A4.R")


input.edits.dt <- fread("result/S71_5__filter_for_A_to_G_sites_for_control/210203-GSE144296.A375-RNA-with-DNA-37-37/210203-GSE144296.A375-RNA-with-DNA-37-37/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")


GSE144296.paired.sequenced.info.dt <- {
    ## read GSE144296 metadat
    fread("external/NCBI.SRA.MetaData/GSE144296.txt") -> .;
    ##s get A375 records
    .[Cell_Line=="A375"] ->.;
    ##s get DNA-RNA-paired samples
    .[, cell_ID_occurrence:=.N, list(cell_ID)][cell_ID_occurrence==2] ->.;
    ##dcast 
    dcast(data=., formula=cell_ID ~ LibrarySource, value.var="Sample Name") ->.;
    ##s select samples with RNA edits identified
    .[TRANSCRIPTOMIC %in% input.edits.dt[, SAMPLE]] ->.;
}


all.DNA.alignments.dt <- foreach(temp.row.dt=iter(GSE144296.paired.sequenced.info.dt, by="row")) %do% {
    cat(date(), "reading DNA alignment for sample ", temp.row.dt[1, GENOMIC], "\n")
    temp.DNA.alignment.bcf.filename <- paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-37-37/210203-GSE144296.A375-DNA-37-37/", temp.row.dt[1, GENOMIC], "/__merged__/DNTRSeq-DNA-trimming/hg38.fa/32/bwa-index-default/DNA/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.bcf")
    temp.DNA.alignment.dt <- fread(cmd=paste(sep="", "bcftools view ", temp.DNA.alignment.bcf.filename, " | grep -v '^##'"), select=1:7) %>% setnames("#CHROM", "CHROM")
    data.table(temp.DNA.alignment.dt, temp.row.dt)
} %>% rbindlist

fwrite(all.DNA.alignments.dt, "./report.ver2/pipeline.validation/210203-GSE144296.A375/all.DNA.alignments.dt.csv.gz")


filtered.RNA.DNA.comparison.dt <- {
    merge(
        x=input.edits.dt,
        y=all.DNA.alignments.dt[, list(CHROM, POS, SAMPLE=TRANSCRIPTOMIC, is.DNA.detected=TRUE)],
        by=c("CHROM", "POS", "SAMPLE"),
        all.x=TRUE, all.y=FALSE) ->.;
    .[is.na(is.DNA.detected)==TRUE, is.DNA.detected:=FALSE]
}

fwrite(filtered.RNA.DNA.comparison.dt, "./report.ver2/pipeline.validation/210203-GSE144296.A375/filtered.RNA.DNA.comparison.dt.csv.gz")

all.raw.RNA.alignments.dt <- foreach(temp.row.dt=iter(GSE144296.paired.sequenced.info.dt, by="row")) %do% {
    cat(date(), "reading RNA alignment for sample ", temp.row.dt[1, GENOMIC], "\n")
    temp.RNA.alignment.bcf.filename <- paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-37-37/210203-GSE144296.A375-RNA-37-37/", temp.row.dt[1, TRANSCRIPTOMIC], "/__merged__/DNTRSeq-RNA-trimming/hg38.fa/32/bwa-index-10.1038_nmeth.2330/32/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.bcf")
    temp.RNA.alignment.dt <- fread(cmd=paste(sep="", "bcftools view ", temp.RNA.alignment.bcf.filename, " | grep -v '^##'"), select=1:7) %>% setnames("#CHROM", "CHROM")
    data.table(temp.RNA.alignment.dt, temp.row.dt)
} %>% rbindlist

raw.RNA.DNA.comparison.dt <- {
    merge(
        x=all.raw.RNA.alignments.dt[, list(CHROM, POS, SAMPLE=TRANSCRIPTOMIC)],
        y=all.DNA.alignments.dt[, list(CHROM, POS, SAMPLE=TRANSCRIPTOMIC, is.DNA.detected=TRUE)],
        by=c("CHROM", "POS", "SAMPLE"),
        all.x=TRUE, all.y=FALSE) ->.;
    .[is.na(is.DNA.detected)==TRUE, is.DNA.detected:=FALSE]
}

fwrite(raw.RNA.DNA.comparison.dt, "./report.ver2/pipeline.validation/210203-GSE144296.A375/raw.RNA.DNA.comparison.dt.csv.gz")



combined.RNA.DNA.comparison.dt <- rbindlist(list(
        data.table(raw.RNA.DNA.comparison.dt[, list(CHROM, POS, SAMPLE, is.DNA.detected)], RNA.edit.type="raw"),
        data.table(filtered.RNA.DNA.comparison.dt[, list(CHROM, POS, SAMPLE, is.DNA.detected)], RNA.edit.type="filtered")))


fwrite(combined.RNA.DNA.comparison.dt, "./report.ver2/pipeline.validation/210203-GSE144296.A375/combined.RNA.DNA.comparison.dt.csv.gz")


{
    
    combined.RNA.DNA.comparison.dt ->.;
    ## get counts of variants per sample
    .[, list(count=.N), list(SAMPLE, is.DNA.detected, RNA.edit.type)] ->.;
    ## add all missing combinations
    setkey(., SAMPLE, is.DNA.detected, RNA.edit.type)
    .[CJ(SAMPLE, is.DNA.detected, RNA.edit.type, unique=TRUE)] ->.;
    setnafill(x=., type="const", fill=0, cols="count")
    ## prettify is.DNA.detected
    .[, DNA.support:=c("Not overlapping with genomic variants", "Overlapping with genomic variants")[(is.DNA.detected == TRUE) + 1]] -> .
    ## take the mean across samples
    .[, list(mean.count=mean(count)), list(DNA.support, RNA.edit.type)] -> .;
    ##
    . -> ..to.plot.dt
    
    ## start plotting
    ggplot(..to.plot.dt, aes(x="a", y=mean.count, fill=DNA.support)) ->.;
    ## plot boxplots
    . + geom_bar(stat="identity", position="dodge") ->.;
    ## facet
    . + facet_grid(RNA.edit.type~., scales="free") ->.;
    ## add theme
    . + theme_pubr(base_size=10) ->.;
    . + theme(axis.text.x=element_blank()) ->.;
    . + guides(fill=guide_legend(nrow=2)) -> .;
    . + labs(x="A>G and A>G;T>C variants", y="Mean(# variants) across samples", fill="") ->.;
    ## save image
    ggsave.A4(
        filename="./report.ver2/pipeline.validation/210203-GSE144296.A375/validation.resized.bar.plot.png",
        plot=.,
        width.r=0.35, height.r=0.25)
    
}
