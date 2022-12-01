library("data.table")
library("readxl")
library("glue")
source("./scripts/common/ggpubr.A4.R")


input.variants.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz")


## A-to-G ratio across simple nucleotide changes
{
    
    input.variants.dt -> .;
    ##c get variant count per sample x variant x Alu type
    .[, list(count=.N), list(SAMPLE, stage, event.summary)] ->.;
    ##s keep simple nucleotide changes only
    .[event.summary %in% c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G")] -> .;
    ##= get variant ratio per sample
    .[, total.count:=sum(count), list(SAMPLE, stage)]
    .[, percentage:=count/total.count]
    .[, percentage.to.plot:=percentage*100]
    ##s get A-to-G variants only
    .[event.summary %in% c('A>G')] ->.;
    
    ##c get A-to-G ratio
    .[, list(A.to.G.ratio=percentage.to.plot), list(SAMPLE, stage)] ->.;
    ## add category
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))] ->.;
    ##
    . -> ..to.plot.dt
    ##
    fwrite(..to.plot.dt[, list(SAMPLE, "A-to-G proportion (%)"=A.to.G.ratio)], "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/all.edits.A.to.G.ratio.across.simple.nucleotide.changes.csv")
    
    ##plot A-to-G ratio
    ggplot(..to.plot.dt, aes(x="a", y=A.to.G.ratio, color=category.ordered)) ->.;
    . + geom_boxplot() ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right") ->.;
    . + guides(color=guide_legend(nrow=3)) ->.;
    ##. + scale_y_continuous(limits=c(0, 100)) ->.;
    . + labs(x="", y="% A-to-G variants", color="") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/all.edits.A.to.G.ratio.across.simple.nucleotide.changes.png",
        plot=.,
        width.r=0.45, height.r=0.15
    )
    
}


