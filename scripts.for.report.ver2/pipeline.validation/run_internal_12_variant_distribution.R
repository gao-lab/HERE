library("data.table")
library("ggpubr")
library("foreach")
library("iterators")
library("glue")
library("stringr")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/"
dir.create(output.directory, recursive=TRUE)

all.variants.dt <- fread("./result/S71_5__filter_for_A_to_G_sites_for_control/210203-GSE144296.A375-RNA-with-DNA-37-37/210203-GSE144296.A375-RNA-with-DNA-37-37/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz")


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
    .[TRANSCRIPTOMIC %in% all.variants.dt[, SAMPLE]] ->.;
}


all.DNA.alignments.dt <- foreach(temp.row.dt=iter(GSE144296.paired.sequenced.info.dt, by="row")) %do% {
    cat(date(), "reading DNA alignment for sample ", temp.row.dt[1, GENOMIC], "\n")
    temp.DNA.alignment.bcf.filename <- paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-37-37/210203-GSE144296.A375-DNA-37-37/", temp.row.dt[1, GENOMIC], "/__merged__/DNTRSeq-DNA-trimming/hg38.fa/32/bwa-index-default/DNA/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.bcf")
    temp.DNA.alignment.dt <- fread(cmd=paste(sep="", "bcftools view ", temp.DNA.alignment.bcf.filename, " | grep -v '^##'"), select=1:7) %>% setnames("#CHROM", "CHROM")
    data.table(temp.DNA.alignment.dt, temp.row.dt)
} %>% rbindlist

fwrite(all.DNA.alignments.dt, glue("{output.directory}/all.DNA.alignments.dt.csv.gz"))


filtered.RNA.DNA.comparison.dt <- {
    merge(
        x=all.variants.dt,
        y=all.DNA.alignments.dt[, list(CHROM, POS, SAMPLE=TRANSCRIPTOMIC, is.DNA.detected=TRUE)],
        by=c("CHROM", "POS", "SAMPLE"),
        all.x=TRUE, all.y=FALSE) ->.;
    .[is.na(is.DNA.detected)==TRUE, is.DNA.detected:=FALSE]
}

fwrite(filtered.RNA.DNA.comparison.dt, glue("{output.directory}/filtered.RNA.DNA.comparison.dt.csv.gz"))

all.raw.RNA.alignments.dt <- foreach(temp.row.dt=iter(GSE144296.paired.sequenced.info.dt, by="row")) %do% {
    cat(date(), "reading RNA alignment for sample ", temp.row.dt[1, GENOMIC], "\n")
    temp.RNA.alignment.bcf.filename <- paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-37-37/210203-GSE144296.A375-RNA-37-37/", temp.row.dt[1, TRANSCRIPTOMIC], "/__merged__/DNTRSeq-RNA-trimming/hg38.fa/32/bwa-index-10.1038_nmeth.2330/32/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.bcf")
    temp.RNA.alignment.dt <- fread(cmd=paste(sep="", "bcftools view ", temp.RNA.alignment.bcf.filename, " | grep -v '^##'"), select=1:7) %>% setnames("#CHROM", "CHROM")
    data.table(temp.RNA.alignment.dt, temp.row.dt)
} %>% rbindlist

## annotate variants of all RNA alignments
temp.all.RNA.alignments.vcf.filename <- glue("{output.directory}/all.RNA.alignments.vcf")
writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##contig=<ID=chr10>",
    "##contig=<ID=chr12>",
    "##contig=<ID=chr14>",
    "##contig=<ID=chr17>",
    "##contig=<ID=chr2>",
    "##contig=<ID=chr3>",
    "##contig=<ID=chr4>",
    "##contig=<ID=chr6>",
    "##contig=<ID=chr7>",
    "##contig=<ID=chr8>",
    "##contig=<ID=chr9>",
    "##contig=<ID=chr11>",
    "##contig=<ID=chr13>",
    "##contig=<ID=chr15>",
    "##contig=<ID=chr16>",
    "##contig=<ID=chr18>",
    "##contig=<ID=chr19>",
    "##contig=<ID=chr20>",
    "##contig=<ID=chr5>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ), temp.all.RNA.alignments.vcf.filename)
fwrite(unique(all.raw.RNA.alignments.dt[, list(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO=".")]), temp.all.RNA.alignments.vcf.filename, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE)
##
temp.all.RNA.alignments.snpEff.vcf.filename <- glue("{output.directory}/all.RNA.alignments.snpEff.vcf")
system(glue("snpEff ann -Xmx30G  -lof -config result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/hg38.fa/32/snpEff.config hg38.fa.GENCODE.32 {temp.all.RNA.alignments.vcf.filename} > {temp.all.RNA.alignments.snpEff.vcf.filename} "))
##
all.RNA.alignments.snpEff.vcf.reformatted.dt <- setnames(fread(cmd=paste(sep="", "cat ", temp.all.RNA.alignments.snpEff.vcf.filename, " | grep -v '^##'"), header=TRUE, drop=c("QUAL", "FILTER")), "#CHROM", "CHROM")[, `:=`(SUBSET=sub(pattern=".*SUBSET=([^;]+);.*", replacement="\\1", x=INFO), ANN=sub(pattern=".*ANN=([^;]+)($|;.*)", replacement="\\1", x=INFO))][grepl("LOF=", INFO) == TRUE, `:=`(LOF=sub(pattern=".*LOF=([^;]+)($|;.*)", replacement="\\1", x=INFO))][, INFO:=NULL]
fwrite(all.RNA.alignments.snpEff.vcf.reformatted.dt, glue("{output.directory}/all.RNA.alignments.snpEff.vcf.reformatted.dt.gz"))
##
all.RNA.alignments.snpEff.vcf.reformatted.ANN.split.dt <- all.RNA.alignments.snpEff.vcf.reformatted.dt[, list(ANN.single.match=str_split(string=ANN, pattern=",")[[1]]), list(CHROM, POS, ID, REF, ALT, ANN)] 
##
ANN.single.match.split.matrix <- do.call(rbind, str_split(string=all.RNA.alignments.snpEff.vcf.reformatted.ANN.split.dt[, ANN.single.match], pattern="\\|"))
all.RNA.alignments.snpEff.vcf.reformatted.ANN.split.annotation.dt <- setnames(data.table(all.RNA.alignments.snpEff.vcf.reformatted.ANN.split.dt[, list(CHROM, POS, ID, REF, ALT)], ANN.single.match.split.matrix), 6:21, c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO"))
all.RNA.alignments.snpEff.vcf.reformatted.ANN.split.annotation.dt[Feature_Type == "transcript", event:=sub(pattern=".*[0-9]+", replacement="\\1", x=HGVS.c)]
all.RNA.alignments.snpEff.vcf.reformatted.ANN.split.annotation.dt[, event:=sub(pattern=".*[0-9]+", replacement="", x=HGVS.c)]

all.RNA.alignments.snpEff.vcf.reformatted.event.summary.dt <- all.RNA.alignments.snpEff.vcf.reformatted.ANN.split.annotation.dt[, list(event.summary=paste(collapse=";", sort(unique(event)))), list(CHROM, POS, ID, REF, ALT)]
fwrite(all.RNA.alignments.snpEff.vcf.reformatted.event.summary.dt, glue("{output.directory}/all.RNA.alignments.snpEff.vcf.reformatted.event.summary.dt.gz"))
##

merge(
    x=all.raw.RNA.alignments.dt[, list(CHROM, POS, REF, ALT, SAMPLE=TRANSCRIPTOMIC)],
    y=all.RNA.alignments.snpEff.vcf.reformatted.event.summary.dt,
    by=c("CHROM", "POS", "REF", "ALT"),
    all.x=TRUE, all.y=FALSE
) -> all.raw.RNA.alignments.with.event.annotation.dt


## build the raw.RNA.DNA.comparison table
raw.RNA.DNA.comparison.dt <- {
    merge(
        x=all.raw.RNA.alignments.with.event.annotation.dt,
        y=all.DNA.alignments.dt[, list(CHROM, POS, SAMPLE=TRANSCRIPTOMIC, is.DNA.detected=TRUE)],
        by=c("CHROM", "POS", "SAMPLE"),
        all.x=TRUE, all.y=FALSE) ->.;
    .[is.na(is.DNA.detected)==TRUE, is.DNA.detected:=FALSE]
}
fwrite(raw.RNA.DNA.comparison.dt, glue("{output.directory}/raw.RNA.DNA.comparison.dt.csv.gz"))

combined.RNA.DNA.comparison.dt <- rbindlist(list(
        data.table(raw.RNA.DNA.comparison.dt[, list(CHROM, POS, SAMPLE, is.DNA.detected, event.summary)], RNA.edit.type="raw"),
        data.table(filtered.RNA.DNA.comparison.dt[, list(CHROM, POS, SAMPLE, is.DNA.detected, event.summary)], RNA.edit.type="filtered")))

fwrite(combined.RNA.DNA.comparison.dt, glue("{output.directory}/combined.RNA.DNA.comparison.dt.csv.gz"))


{
    
    copy(combined.RNA.DNA.comparison.dt) ->.;
    ## discard variants not overlapping with known transcripts
    .[(event.summary %in% c(NA, "")) == FALSE] -> .;    
    ## collapse mixed events (i.e., those that are more than plus/minus strand variants)
    .[, event.summary.collapsed:=c(
            ## single strand
            "A>C"="A>C", "A>G"="A>G", "A>T"="A>T",
            "C>A"="C>A", "C>G"="C>G", "C>T"="C>T",
            "G>A"="G>A", "G>C"="G>C", "G>T"="G>T",
            "T>A"="T>A", "T>C"="T>C", "T>G"="T>G",
            ## plus + minus strand
            "A>C;T>G"="A>C;T>G", "A>G;T>C"="A>G;T>C", "A>T;T>A"="A>T;T>A",
            "C>A;G>T"="C>A;G>T", "C>G;G>C"="C>G;G>C", "C>T;G>A"="C>T;G>A"
        )[event.summary]] -> .;
    .[is.na(event.summary.collapsed) == TRUE, event.summary.collapsed:="mixed variants"] -> .;
    ##
    ## get counts of variants per sample
    .[, list(count=.N), list(SAMPLE, is.DNA.detected, RNA.edit.type, event.summary.collapsed)] ->.;
    ## add all missing combinations
    setkey(., SAMPLE, is.DNA.detected, RNA.edit.type, event.summary.collapsed)
    .[CJ(SAMPLE, is.DNA.detected, RNA.edit.type, event.summary.collapsed, unique=TRUE)] ->.;
    setnafill(x=., type="const", fill=0, cols="count")
    ## prettify is.DNA.detected
    .[, DNA.support:=c("Not overlapping with genomic variants", "Overlapping with genomic variants")[(is.DNA.detected == TRUE) + 1]] -> .
    ## prettify event.summary.collapsed
    .[, event.summary.collapsed.ordered:=factor(
            event.summary.collapsed,
            levels=c(
                "A>C",
                "A>G",
                "A>T",
                "C>A",
                "C>G",
                "C>T",
                "G>A",
                "G>T",
                "G>C",
                "T>A",
                "T>C",
                "T>G",
                "A>C;T>G",
                "A>G;T>C",
                "A>T;T>A",
                "C>A;G>T",
                "C>G;G>C",
                "C>T;G>A",
                "mixed variants"
            ))] -> .;
    ## take the mean across samples
    .[, list(mean.count=mean(count)), list(event.summary.collapsed.ordered, DNA.support, RNA.edit.type)] -> .;
    ## 
    . -> ..to.plot.dt
    
    ## start plotting
    ggplot(..to.plot.dt[(grepl(";", event.summary.collapsed.ordered) == FALSE) & (event.summary.collapsed.ordered != "mixed variants"), ], aes(x=event.summary.collapsed.ordered, y=mean.count, fill=DNA.support)) ->.;
    ## plot boxplots
    . + geom_bar(stat="identity", position="dodge") ->.;
    ## facet
    . + facet_grid(RNA.edit.type~., scales="free_y") ->.;
    ## add theme
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) ->.;
    . + labs(x="Mismatch type", y="Mean(# variants) across samples", fill="") ->.;
    ## save image
    ggsave.A4(
        filename=glue("{output.directory}/validation.all.variants.unambiguous.strand.only.png"),
        plot=.,
        width.r=0.9, height.r=0.45)

    ggplot(..to.plot.dt, aes(x=event.summary.collapsed.ordered, y=mean.count, fill=DNA.support)) ->.;
    ## plot boxplots
    . + geom_bar(stat="identity", position="dodge") ->.;
    ## facet
    . + facet_grid(RNA.edit.type~., scales="free_y") ->.;
    ## add theme
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) ->.;
    . + labs(x="Mismatch type", y="Mean(# variants) across samples", fill="") ->.;
    ## save image
    ggsave.A4(
        filename=glue("{output.directory}/validation.all.variants.all.png"),
        plot=.,
        width.r=0.9, height.r=0.45)
    
}

{
    
    copy(combined.RNA.DNA.comparison.dt) ->.;
    ## discard variants not overlapping with known transcripts
    .[(event.summary %in% c(NA, "")) == FALSE] -> .;    
    ## collapse mixed events (i.e., those that are more than plus/minus strand variants)
    .[, event.summary.collapsed:=c(
            ## single strand
            "A>C"="A>C", "A>G"="A>G", "A>T"="A>T",
            "C>A"="C>A", "C>G"="C>G", "C>T"="C>T",
            "G>A"="G>A", "G>C"="G>C", "G>T"="G>T",
            "T>A"="T>A", "T>C"="T>C", "T>G"="T>G",
            ## plus + minus strand
            "A>C;T>G"="A>C;T>G", "A>G;T>C"="A>G;T>C", "A>T;T>A"="A>T;T>A",
            "C>A;G>T"="C>A;G>T", "C>G;G>C"="C>G;G>C", "C>T;G>A"="C>T;G>A"
        )[event.summary]] -> .;
    .[is.na(event.summary.collapsed) == TRUE, event.summary.collapsed:="mixed variants"] -> .;
    ##
    ## get counts of variants per sample
    .[, list(count=.N), list(SAMPLE, is.DNA.detected, RNA.edit.type, event.summary.collapsed)] ->.;
    ## add all missing combinations
    setkey(., SAMPLE, is.DNA.detected, RNA.edit.type, event.summary.collapsed)
    .[CJ(SAMPLE, is.DNA.detected, RNA.edit.type, event.summary.collapsed, unique=TRUE)] ->.;
    setnafill(x=., type="const", fill=0, cols="count")
    ## prettify is.DNA.detected
    .[, DNA.support:=c("Not overlapping with genomic variants", "Overlapping with genomic variants")[(is.DNA.detected == TRUE) + 1]] -> .
    ## prettify event.summary.collapsed
    .[, event.summary.collapsed.ordered:=factor(
            event.summary.collapsed,
            levels=c(
                "A>C",
                "A>G",
                "A>T",
                "C>A",
                "C>G",
                "C>T",
                "G>A",
                "G>C",
                "G>T",
                "T>A",
                "T>C",
                "T>G",
                "A>C;T>G",
                "A>G;T>C",
                "A>T;T>A",
                "C>A;G>T",
                "C>G;G>C",
                "C>T;G>A",
                "mixed variants"
            ))] -> .;
    ## take only filtered, DNA undetected edits
    .[RNA.edit.type=='filtered'] -> .;
    .[is.DNA.detected==FALSE] -> .;
    ##
    .[event.summary.collapsed %in% c(
                "A>C",
                "A>G",
                "A>T",
                "C>A",
                "C>G",
                "C>T",
                "G>A",
                "G>T",
                "G>C",
                "T>A",
                "T>C",
                "T>G"
                )] ->.;
    ## take the A-to-G + A-to-G/T-to-C ratio per sample
    .[, total.count:=sum(count), list(SAMPLE)] -> .;
    .[, percentage:=count/total.count] -> .;
    .[, percentage.to.plot:=percentage*100] -> .;
    ## filter for samples with total count >= 10
    .[total.count >= 10] -> .;
    . -> to.plot.raw.dt
    ## compute mean and std
    .[, list(percentage.to.plot.mean=mean(percentage.to.plot), percentage.to.plot.stdev=sqrt(var(percentage.to.plot))), list(event.summary.collapsed.ordered)] -> .;
    . -> ..to.plot.dt

    copy(..to.plot.dt) -> .;
    ggplot(., aes(x=event.summary.collapsed.ordered, y=percentage.to.plot.mean)) ->.;
    . + geom_bar(stat="identity", color="black", fill="white") ->.; # leave white space for point plotting
    . + geom_errorbar(aes(ymin=percentage.to.plot.mean, ymax=percentage.to.plot.mean + percentage.to.plot.stdev)) -> .;
    . + geom_point(data=to.plot.raw.dt, aes(x=event.summary.collapsed.ordered, y=percentage.to.plot)) -> .; # add individual points
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) -> .;
    . + scale_y_continuous(limits=c(0, 100)) ->.;
    . + labs(x="", y="%", color="") ->.;
    ggsave.A4(
        filename=glue("{output.directory}/validation.all.unambiguous.variants.percentage.bar.png"),
        plot=.,
        width.r=0.8, height.r=0.45)

}
