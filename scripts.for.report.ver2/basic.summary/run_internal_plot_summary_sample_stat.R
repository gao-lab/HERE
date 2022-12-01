library("data.table")
library("glue")
source("./scripts/common/ggpubr.A4.R")


combined.summary.dt <- fread("./result/B92_1__summarize_sample_stat/210215-sixth-dataset/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/combined.summary.dt.gz")

fwrite(combined.summary.dt[, list(SAMPLE, "mapping rate"=mapping.rate)], "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/mapping.rates.csv")
fwrite(combined.summary.dt[, list(SAMPLE, "mean depth across whole genome"=mean.depth.across.whole.genome)], "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/mean.depth.across.whole.genome.csv")

## A-to-G ratio across simple nucleotide changes
{
    
    combined.summary.dt -> .;
    . -> ..to.plot.dt
    ##

    ##plot mapping.rate
    ggplot(..to.plot.dt, aes(x="a", y=mapping.rate)) -> .;
    . + geom_boxplot() ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) ->.;
    . + labs(x="", y="Mapping rate", color="") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/mapping.rate.png",
        plot=.,
        width.r=0.25, height.r=0.25
    )
    
    ##plot mean.depth.across.whole.genome
    ggplot(..to.plot.dt, aes(x="a", y=mean.depth.across.whole.genome)) -> .;
    . + geom_boxplot() ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) ->.;
    . + labs(x="", y="Mean depth\nacross whole genome", color="") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/mean.depth.across.whole.genome.png",
        plot=.,
        width.r=0.25, height.r=0.25
    )
    
}


