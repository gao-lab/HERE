library("data.table")
library("readxl")
library("glue")
source("./scripts/common/ggpubr.A4.R")

input.edits.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")

input.variants.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz")


## A-to-G ratio per Alu-or-not, boxplot
{
    input.variants.dt ->.
    ##c get variant count per sample x variant x Alu type
    .[, list(count=.N), list(SAMPLE, stage, event.summary, Alu.or.not=c("non-Alu", "Alu")[(SUBSET=='Alu') + 1])] ->.;
    ##= get variant ratio per sample
    .[, total.count:=sum(count), list(SAMPLE, stage, Alu.or.not)]
    .[, percentage:=count/total.count]
    .[, percentage.to.plot:=percentage*100]
    ##s get A-to-G variants only
    .[event.summary %in% c('A>G', 'A>G;T>C')] ->.;
    ##c get A-to-G ratio
    .[, list(A.to.G.ratio=sum(percentage.to.plot)), list(SAMPLE, stage, Alu.or.not)] ->.;
    ## add category
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))] ->.;
    ##
    . -> ..to.plot.dt
    ##plot A-to-G ratio
    ggplot(..to.plot.dt, aes(x="a", y=A.to.G.ratio, color=category.ordered)) ->.;
    . + geom_boxplot() ->.;
    . + facet_grid(~Alu.or.not) -> .;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right") ->.;
    . + guides(color=guide_legend(nrow=3)) ->.;
    ##. + scale_y_continuous(limits=c(0, 100)) ->.;
    . + labs(x="", y="% A-to-G variants", color="") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/all.edits.A.to.G.ratio.per.Alu.or.not.boxplot.png",
        plot=.,
        width.r=0.45, height.r=0.15
    )       
}


