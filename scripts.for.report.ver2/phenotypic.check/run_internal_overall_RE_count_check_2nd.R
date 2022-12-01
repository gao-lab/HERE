library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("scales")
library("ggdendro")
library("glue")
source("./scripts/common/ggpubr.A4.R")

## read total edits in all samples
edit.info.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")

## read other tables (for filtering)
subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")



normal.RE.and.gene.per.stage.dt <- {
    subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt ->.;
    ##s keep matched stages only
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell')] ->.;
    ##= keep columns needed 
    .[, list(CHROM, POS, Gene_ID, Gene_Name, Annotation.corrected, stage)] %>% unique ->.;
}





edit.normal.RE.and.gene.only.dt <- {
    edit.info.dt ->.;
    ##s keep matched stages only
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', 'zygote.2PN', '2-cell', '4-cell', '8-cell', 'morula')] ->.;
    .[, stage.corrected:=stage][stage=="zygote.2PN", stage.corrected:="zygote"] ->.;
    ##m keep normal-RE edits only
    merge(x=., y=unique(normal.RE.and.gene.per.stage.dt[, list(CHROM, POS, stage)]),
          by.x=c("CHROM", "POS", "stage.corrected"), by.y=c("CHROM", "POS", "stage"),
          all.x=FALSE, all.y=FALSE, allow.cartesian=TRUE) ->.;
}

## compute per sample: count of REs per gene
normal.RE.and.gene.count.info.dt <- edit.normal.RE.and.gene.only.dt[, list(count.RE=length(unique(ID))), list(SAMPLE, gse, stage, stage.corrected, is.normal, treatment, disease, maternal.age)]

{
    
    ## select those with abnormal or old embryos
    normal.RE.and.gene.count.info.dt[gse %in% c("GSE133854", "GSE95477")] ->.;
    ## prettify group labels
    ##-# 1. GSE133854
    .[gse=="GSE133854" & disease == "androgenetic", final.label:="AG"]
    .[gse=="GSE133854" & disease == "parthenogenetic", final.label:="PG"]
    .[gse=="GSE133854" & is.normal == TRUE & stage != "oocyte.MII", final.label:="BI"]
    ##-# 2. GSE95477
    .[gse=="GSE95477" & as.integer(maternal.age) > 35, final.label:="old"]
    .[gse=="GSE95477" & as.integer(maternal.age) < 35, final.label:="young"]
    ##-# final cleaning
    .[is.na(final.label)==FALSE] ->.;
    ## prettify final.label
    .[, final.label.reordered:=factor(final.label, levels=c("old", "young", "amanitin", "control", "AG", "PG", "BI"))]
    ## prettify stages
    temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    . -> to.plot.dt
    ##
    ## start plotting
    common.option <- function(temp.x){temp.x + facet_wrap(~stage.description.ordered, scales="free", nrow=1) + theme_pubr(base_size=11) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous(expand=expansion(mult=c(0.1, 0.2)))}
    (ggplot(.[gse=="GSE95477"], aes(x=final.label.reordered, y=count.RE)) + geom_boxplot() + stat_compare_means(comparisons = list(c("old", "young")), method="wilcox.test", method.args=list(alternative="less"))  + ggtitle("GSE95477") + labs(x="", y="# REE-matching edits")) %>% common.option -> temp.GSE95477.ggplot;
    (ggplot(.[gse=="GSE133854"], aes(x=final.label.reordered, y=count.RE)) + geom_boxplot() + stat_compare_means(comparisons = list(c("PG", "BI"), c("AG", "BI")), method="wilcox.test", method.args=list(alternative="less")) + ggtitle("GSE133854") + labs(x="", y="")) %>% common.option -> temp.GSE133854.ggplot;
    ggarrange(
        plotlist=list(
            temp.GSE95477.ggplot,
            temp.GSE133854.ggplot
        ),
        nrow=1, widths=c(1, 1.75), align="h") ->.;
    ## save image
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/REE.normal.valid.gene.only.REE.count.histogram.three.datasets.png",
        plot=.,
        width.r=0.9, height.r=0.3)

}


## Compute sample size, effect size, and 95CI

{

    copy(to.plot.dt) -> A;
    A[gse=="GSE95477"][stage=='oocyte.GV'] -> A;
    wilcox.test(A[final.label=="old", count.RE], A[final.label=="young", count.RE], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(gse="GSE95477", stage="oocyte (GV)", left.group="old", right.group="young", count.of.samples.for.the.left.group=A[final.label=="old", .N], count.of.samples.for.the.right.group=A[final.label=="young", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.GSE95477.oocyte.GV.dt   

    copy(to.plot.dt) -> A;
    A[gse=="GSE95477"][stage=='oocyte.MII'] -> A;
    wilcox.test(A[final.label=="old", count.RE], A[final.label=="young", count.RE], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(gse="GSE95477", stage="oocyte (MII)", left.group="old", right.group="young", count.of.samples.for.the.left.group=A[final.label=="old", .N], count.of.samples.for.the.right.group=A[final.label=="young", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.GSE95477.oocyte.MII.dt   

    copy(to.plot.dt) -> A;
    A[gse=="GSE133854"][stage=='zygote'] -> A;
    wilcox.test(A[final.label=="AG", count.RE], A[final.label=="BI", count.RE], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(gse="GSE133854", stage="zygote", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[final.label=="AG", .N], count.of.samples.for.the.right.group=A[final.label=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.GSE133854.zygote.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[gse=="GSE133854"][stage=='zygote'] -> A;
    wilcox.test(A[final.label=="PG", count.RE], A[final.label=="BI", count.RE], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(gse="GSE133854", stage="zygote", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[final.label=="PG", .N], count.of.samples.for.the.right.group=A[final.label=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.GSE133854.zygote.PG.vs.BI.dt   


    copy(to.plot.dt) -> A;
    A[gse=="GSE133854"][stage=='2-cell'] -> A;
    wilcox.test(A[final.label=="AG", count.RE], A[final.label=="BI", count.RE], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(gse="GSE133854", stage="2-cell", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[final.label=="AG", .N], count.of.samples.for.the.right.group=A[final.label=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.GSE133854.2cell.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[gse=="GSE133854"][stage=='2-cell'] -> A;
    wilcox.test(A[final.label=="PG", count.RE], A[final.label=="BI", count.RE], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(gse="GSE133854", stage="2-cell", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[final.label=="PG", .N], count.of.samples.for.the.right.group=A[final.label=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.GSE133854.2cell.PG.vs.BI.dt   



    copy(to.plot.dt) -> A;
    A[gse=="GSE133854"][stage=='4-cell'] -> A;
    wilcox.test(A[final.label=="AG", count.RE], A[final.label=="BI", count.RE], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(gse="GSE133854", stage="4-cell", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[final.label=="AG", .N], count.of.samples.for.the.right.group=A[final.label=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.GSE133854.4cell.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[gse=="GSE133854"][stage=='4-cell'] -> A;
    wilcox.test(A[final.label=="PG", count.RE], A[final.label=="BI", count.RE], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(gse="GSE133854", stage="4-cell", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[final.label=="PG", .N], count.of.samples.for.the.right.group=A[final.label=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.GSE133854.4cell.PG.vs.BI.dt   


    copy(to.plot.dt) -> A;
    A[gse=="GSE133854"][stage=='8-cell'] -> A;
    wilcox.test(A[final.label=="AG", count.RE], A[final.label=="BI", count.RE], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(gse="GSE133854", stage="8-cell", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[final.label=="AG", .N], count.of.samples.for.the.right.group=A[final.label=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.GSE133854.8cell.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[gse=="GSE133854"][stage=='8-cell'] -> A;
    wilcox.test(A[final.label=="PG", count.RE], A[final.label=="BI", count.RE], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(gse="GSE133854", stage="8-cell", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[final.label=="PG", .N], count.of.samples.for.the.right.group=A[final.label=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.GSE133854.8cell.PG.vs.BI.dt   


    list(
        temp.GSE95477.oocyte.GV.dt,
        temp.GSE95477.oocyte.MII.dt,
        temp.GSE133854.zygote.AG.vs.BI.dt, temp.GSE133854.zygote.PG.vs.BI.dt,
        temp.GSE133854.2cell.AG.vs.BI.dt, temp.GSE133854.2cell.PG.vs.BI.dt,
        temp.GSE133854.4cell.AG.vs.BI.dt, temp.GSE133854.4cell.PG.vs.BI.dt,
        temp.GSE133854.8cell.AG.vs.BI.dt, temp.GSE133854.8cell.PG.vs.BI.dt
    ) -> .;
    rbindlist(., use.names=TRUE) -> .;

    fwrite(., "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/stat.for.REE.normal.valid.gene.only.REE.count.histogram.three.datasets.csv.gz")
}
