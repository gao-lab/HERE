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

    ## plot editing level per group, histogram, log10-scale
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

    ## plot editing level per group, histogram, no editing level 1, as-is
    ggplot(..to.plot.dt[AF != 1], aes(x=AF)) ->.;
    . + geom_histogram() ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) ->.;
    . + facet_wrap(~stage.description.ordered, scales="free_y", ncol=3) -> .;
    . + labs(x="Editing level of all those edits whose editing level is not 1\nacross all samples of the same group\n(lower-bounded by 0.1 due to identification pipeline)", y="Count of such (edit, sample)") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/all.editing.levels.per.stage.no.AF.1.histogram.yaxis.asis.png",
        plot=.,
        width.r=0.9, height.r=0.9
    )    
    
}
