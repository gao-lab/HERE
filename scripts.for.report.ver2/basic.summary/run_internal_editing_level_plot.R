library("data.table")
library("readxl")
library("glue")
source("./scripts/common/ggpubr.A4.R")

input.edits.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")


## per-stage count, barplot
{
    
    input.edits.dt ->.
    ## prettify stage
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    . -> ..to.plot.dt

    ## plot editing level per group, histogram
    ggplot(..to.plot.dt, aes(x=AF)) ->.;
    . + geom_histogram() ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) ->.;
    . + facet_wrap(~stage.description.ordered, scales="free_y", ncol=3) -> .;
    . + scale_y_log10( labels=function(x) format(x, scientific=FALSE)) -> .;
    . + labs(x="Editing level of all edits across all samples of the same group\n(lower-bounded by 0.1 due to identification pipeline)", y="Count of such (edit, sample)\n(in log10 scale)") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/all.editing.levels.per.stage.histogram.png",
        plot=.,
        width.r=0.9, height.r=0.9
    )

    
    ## deprecated: too messy, will try per-facet histogram instead
    ## plot editing level per group, boxplot
    ggplot(..to.plot.dt, aes(x=stage.description.ordered, y=AF)) ->.;
    . + geom_boxplot() ->.;
    . + coord_flip() ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) ->.;
    . + scale_y_continuous(limits=c(0.1, 1), breaks=seq(0.1, 1, 0.3)) -> .;
    . + labs(x="", y="editing level\n(lower-bounded by 0.1 due to identification pipeline)") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/all.editing.levels.per.stage.barplot.png",
        plot=.,
        width.r=0.9, height.r=0.7
    )

    ## violin style, hardly any signal
    ggplot(..to.plot.dt, aes(x=stage.description.ordered, y=AF)) ->.;
    . + geom_violin() ->.;
    . + coord_flip() ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) ->.;
    . + scale_y_continuous(limits=c(0.1, 1), breaks=seq(0.1, 1, 0.3)) -> .;
    . + labs(x="", y="editing level\n(lower-bounded by 0.1 due to identification pipeline)") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/all.editing.levels.per.stage.violin.png",
        plot=.,
        width.r=0.9, height.r=0.7
    )
    
}
