library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("scales")
library("ggdendro")
source("./scripts/common/ggpubr.A4.R")

## read total edits in all samples
edit.info.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")

## read other tables (for filtering)
## 1. REs
## subset.recurrent.edits.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.dt.txt.gz")
## 2. RE and their genes
subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

## normal.RE.per.stage.dt <- {
##     subset.recurrent.edits.only.dt ->.;
##     ##s keep matched stages only
##     .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula')] ->.;
##     ##= keep columns needed 
##     .[, list(CHROM, POS, stage)] %>% unique ->.;
## }


normal.RE.and.gene.per.stage.dt <- {
    subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt ->.;
    ##s keep matched stages only
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell')] ->.;
    ##= keep columns needed 
    .[, list(CHROM, POS, Gene_ID, Gene_Name, Annotation.corrected, stage)] %>% unique ->.;
}


## edit.normal.RE.only.dt <- {
##     edit.info.dt ->.;
##     ##s keep matched stages only
##     .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', 'zygote.2PN', '2-cell', '4-cell', '8-cell', 'morula')] ->.;
##     .[, stage.corrected:=stage][stage=="zygote.2PN", stage.corrected:="zygote"] ->.;
##     ##m keep normal-RE edits only
##     merge(x=., y=normal.RE.per.stage.dt,
##           by.x=c("CHROM", "POS", "stage.corrected"), by.y=c("CHROM", "POS", "stage"),
##           all.x=FALSE, all.y=FALSE, allow.cartesian=TRUE) ->.;
## }

## ## compute per sample: count of REs per gene
## normal.RE.count.info.dt <- edit.normal.RE.only.dt[, list(count.RE=length(unique(ID))), list(SAMPLE, gse, stage, stage.corrected, is.normal, treatment, disease, maternal.age)]

## {
##     ## select those with abnormal or old embryos
##     normal.RE.count.info.dt[gse %in% c("GSE101571", "GSE133854", "GSE95477", "GSE65481")] ->.;
##     ## prettify stages
##     temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
##     merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
##     ## start plotting
##     ggplot(., aes(x=paste(sep="", disease, ";", treatment, "", maternal.age > 35), y=count.RE)) ->.;
##     ## plot boxplot
##     . + geom_boxplot() ->.;
##     ## facet
##     . + facet_wrap(~gse + stage, scales="free") ->.;
##     ## add theme
##     . + theme_pubr() ->.;
##     . + theme(axis.text.x = element_text(angle=45, hjust=1)) ->.;
## ##    . + scale_y_continuous(limits=c(0, 6.5)) ->.;
##     . + labs(x="", y="# RE-matching edits identified") ->.;
## ##    . + ggtitle("All 2-cells") ->.;
##     ## save image
##     ggsave.A4(
##         filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/normal.RE.count.histogram.four.datasets.png",
##         plot=.,
##         width.r=3, height.r=1.5)
## }








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
    normal.RE.and.gene.count.info.dt[gse %in% c("GSE101571", "GSE133854", "GSE95477", "GSE65481")] ->.;
    ## prettify group labels
    ##-# 1. GSE101571
    .[gse=="GSE101571" & treatment == "control", final.label:="control"]
    .[gse=="GSE101571" & treatment == "amanitin", final.label:="amanitin"]
    ##-# 2. GSE133854
    .[gse=="GSE133854" & disease == "androgenetic", final.label:="AG"]
    .[gse=="GSE133854" & disease == "parthenogenetic", final.label:="PG"]
    .[gse=="GSE133854" & is.normal == TRUE & stage != "oocyte.MII", final.label:="BI"]
    ##-# 3. GSE65481
    .[gse=="GSE65481" & disease == "viability predicted to be bad", final.label:="bad"]
    .[gse=="GSE65481" & disease == "viability predicted to be good", final.label:="good"]
    ##-# 4. GSE95477
    .[gse=="GSE95477" & as.integer(maternal.age) > 35, final.label:="old"]
    .[gse=="GSE95477" & as.integer(maternal.age) < 35, final.label:="young"]
    ##-# final cleaning
    .[is.na(final.label)==FALSE] ->.;
    ## prettify final.label
    .[, final.label.reordered:=factor(final.label, levels=c("old", "young", "bad", "good", "amanitin", "control", "AG", "PG", "BI"))]
    ## prettify stages
    temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## start plotting
    common.option <- function(temp.x){temp.x + facet_wrap(~stage.description.ordered, scales="free", nrow=1) + theme_pubr() + theme(axis.text.x = element_text(angle=45, hjust=1)) + scale_y_continuous(expand=expansion(mult=c(0.1, 0.2)))}
    (ggplot(.[gse=="GSE95477"], aes(x=final.label.reordered, y=count.RE)) + geom_boxplot() + stat_compare_means(comparisons = list(c("old", "young")), method="wilcox.test", method.args=list(alternative="less"))  + ggtitle("GSE95477") + labs(x="", y="# RE-matching edits identified\non protein-coding genes")) %>% common.option -> temp.GSE95477.ggplot;
    (ggplot(.[gse=="GSE65481"], aes(x=final.label.reordered, y=count.RE)) + geom_boxplot() + stat_compare_means(comparisons = list(c("bad", "good")), method="wilcox.test", method.args=list(alternative="less")) + ggtitle("GSE65481") + labs(x="", y="")) %>% common.option -> temp.GSE65481.ggplot;
    (ggplot(.[gse=="GSE133854"], aes(x=final.label.reordered, y=count.RE)) + geom_boxplot() + stat_compare_means(comparisons = list(c("AG", "BI"), c("PG", "BI")), method="wilcox.test", method.args=list(alternative="less")) + ggtitle("GSE133854") + labs(x="", y="# RE-matching edits identified\non protein-coding genes")) %>% common.option -> temp.GSE133854.ggplot;
    (ggplot(.[gse=="GSE101571"], aes(x=final.label.reordered, y=count.RE)) + geom_boxplot() + stat_compare_means(comparisons = list(c("amanitin", "control")), method="wilcox.test", method.args=list(alternative="less")) + ggtitle("GSE101571") + labs(x="", y="")) %>% common.option -> temp.GSE101571.ggplot;
    ggarrange(
        plotlist=list(
            ggarrange(plotlist=list(temp.GSE95477.ggplot, temp.GSE65481.ggplot, temp.GSE101571.ggplot), nrow=1, widths=c(4,2,3), align="h"),
            temp.GSE133854.ggplot),
        nrow=2) ->.;
    ##. + labs(x="", y="# RE-matching edits identified\non protein-coding genes") ->.;
##    . + ggtitle("All 2-cells") ->.;
    ## save image
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/normal.valid.gene.only.RE.count.histogram.four.datasets.png",
        plot=.,
        width.r=0.9, height.r=0.5)
}
