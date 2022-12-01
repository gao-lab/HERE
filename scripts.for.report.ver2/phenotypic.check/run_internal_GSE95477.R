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
## 1. RE and their genes
subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")
## 2. annotation of maternal genes
combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt <- fread("result/S42_1__annotate_embryonic_genes/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/201221-fifth-phenotype-collection/combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt.txt.gz")


normal.RE.and.gene.per.stage.dt <- {
    subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt ->.;
    ##s keep matched stages only
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell')] ->.;
    ##= keep columns needed 
    .[, list(CHROM, POS, Gene_ID, Gene_Name, Annotation.corrected, stage)] %>% unique ->.;
}

maternal.gene.info.dt <- combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt[cluster %in% c("maternal.decay", "maternal.others"), list(Gene_ID=gene.id, Gene_Name=gene.name, cluster)]

edit.normal.RE.only.maternal.gene.only.dt <- {
    edit.info.dt ->.;
    ##s keep matched stages only
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', 'zygote.2PN', '2-cell', '4-cell', '8-cell')] ->.;
    .[, stage.corrected:=stage][stage=="zygote.2PN", stage.corrected:="zygote"] ->.;
    ##m keep normal-RE edits only
    merge(x=., y=normal.RE.and.gene.per.stage.dt,
          by.x=c("CHROM", "POS", "stage.corrected"), by.y=c("CHROM", "POS", "stage"),
          all.x=FALSE, all.y=FALSE, allow.cartesian=TRUE) ->.;
    ## keep maternal genes only
    merge(x=., y=maternal.gene.info.dt,
          by=c("Gene_ID", "Gene_Name"),
          all.x=FALSE, all.y=FALSE) ->.;
}

## compute per sample: count of REs per gene
RE.count.info.dt <- edit.normal.RE.only.maternal.gene.only.dt[, list(count.RE=length(unique(ID)), count.gene=length(unique(Gene_ID)), mean.RE.per.gene=length(unique(ID)) / length(unique(Gene_ID)) ), list(SAMPLE, gse, stage, stage.corrected, is.normal, treatment, disease, maternal.age, cluster)]

## Plot all samples of GSE95477
{
    
    RE.count.info.dt[gse=="GSE95477"] ->.;
    ## bin maternal ages
    .[, maternal.bin.ordered:=factor(c("young", "old")[(as.integer(maternal.age) > 35)+1], levels=c("old", "young")) ]
    ##= prettify cluster
    .[, cluster.renamed.ordered:=factor(c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster])]
    ## prettify stages
    temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    . -> REE.count.boxplot.GSE95477.to.plot.dt
    
    ## start plotting
    REE.count.boxplot.GSE95477.to.plot.dt -> .;
    ggplot(., aes(x=maternal.bin.ordered, y=mean.RE.per.gene, fill=maternal.bin.ordered)) ->.;
    ## add mean comparison
    . + stat_compare_means(comparisons=list(c("old", "young")), method="wilcox.test", method.args=list(alternative="less")) ->.;
    ## plot boxplots
    . + geom_boxplot() ->.;
    ## facet
    . + facet_grid(cluster.renamed.ordered~stage.description.ordered, scales="free_y") ->.;
    ## add theme
    . + theme_pubr() ->.;
    . + theme(axis.text.x = element_blank()) ->.;
    . + scale_y_continuous(limits=c(4.5, 6.5)) ->.;
    . + labs(x="", y="Average # REE-matching edits per gene", fill="Maternal age group") ->.;
    ## save image
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/REE.count.boxplot.GSE95477.png",
        plot=.,
        width.r=0.45, height.r=0.4)
}

{

    copy(REE.count.boxplot.GSE95477.to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='oocyte.GV'] -> A;
    wilcox.test(A[maternal.bin.ordered=="old", mean.RE.per.gene], A[maternal.bin.ordered=="young", mean.RE.per.gene], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="decay at 8-cell", stage="oocyte (GV)", left.group="old", right.group="young", count.of.samples.for.the.left.group=A[maternal.bin.ordered=="old", .N], count.of.samples.for.the.right.group=A[maternal.bin.ordered=="young", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.decay.at.8.cell.oocyte.GV.dt   


    copy(REE.count.boxplot.GSE95477.to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='oocyte.GV'] -> A;
    wilcox.test(A[maternal.bin.ordered=="old", mean.RE.per.gene], A[maternal.bin.ordered=="young", mean.RE.per.gene], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="oocyte (GV)", left.group="old", right.group="young", count.of.samples.for.the.left.group=A[maternal.bin.ordered=="old", .N], count.of.samples.for.the.right.group=A[maternal.bin.ordered=="young", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.others.oocyte.GV.dt   


    copy(REE.count.boxplot.GSE95477.to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='oocyte.MII'] -> A;
    wilcox.test(A[maternal.bin.ordered=="old", mean.RE.per.gene], A[maternal.bin.ordered=="young", mean.RE.per.gene], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="decay at 8-cell", stage="oocyte (MII)", left.group="old", right.group="young", count.of.samples.for.the.left.group=A[maternal.bin.ordered=="old", .N], count.of.samples.for.the.right.group=A[maternal.bin.ordered=="young", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.decay.at.8.cell.oocyte.MII.dt   


    copy(REE.count.boxplot.GSE95477.to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='oocyte.MII'] -> A;
    wilcox.test(A[maternal.bin.ordered=="old", mean.RE.per.gene], A[maternal.bin.ordered=="young", mean.RE.per.gene], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="oocyte (MII)", left.group="old", right.group="young", count.of.samples.for.the.left.group=A[maternal.bin.ordered=="old", .N], count.of.samples.for.the.right.group=A[maternal.bin.ordered=="young", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.others.oocyte.MII.dt   

    list(
        temp.decay.at.8.cell.oocyte.GV.dt,
        temp.others.oocyte.GV.dt,
        temp.decay.at.8.cell.oocyte.MII.dt,
        temp.others.oocyte.MII.dt
    ) -> .;
    rbindlist(., use.names=TRUE) -> .;

    fwrite(., "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/stat.for.REE.count.boxplot.GSE95477.png.csv.gz")
}




## Plot clustering result

{
    edit.normal.RE.only.maternal.gene.only.dt ->.;
    ##s select GSE95477 only
    .[gse %in% c("GSE95477")] ->.;
    ##c get RE (on maternal genes) per sample and cluster
    .[, list(SAMPLE, stage, maternal.age, maternal.bin.ordered=factor(c("young", "old")[(as.integer(maternal.age) > 35)+1], levels=c("old", "young")))] %>% unique -> temp.phenotype.dt
    .[, list(CHROM, POS, SAMPLE, one=1)] %>% unique ->.;
    ##d dcast
    dcast(., CHROM + POS ~ SAMPLE, value.var="one", fill=0) ->.;
    ##t transpose
    t(.[, c(-1, -2)]) -> .;
    ## generate cluster plot
    dendro_data(hclust(dist(., method="manhattan"))) ->.;
    .[["labels"]] <- merge(x=.[["labels"]], y=temp.phenotype.dt, by.x="label", by.y="SAMPLE", all=TRUE)
    ggplot() + geom_segment(data=segment(.), aes(x=x, y=y, xend=xend, yend=yend)) + geom_text(data=label(.), aes(x=x, y=y, label=paste(sep="", label, " ; age: ", maternal.age), hjust=0, color=paste(sep="", stage, ";", maternal.bin.ordered))) + coord_flip() +  scale_y_continuous(expand=expansion(mult=c(0.3, 1.3)), breaks=seq(0, 3000, 1000), trans=reverse_trans()) -> .;
    . + theme_pubr() ->.;
    . + theme(axis.text.y=element_blank(), axis.ticks.y = element_blank()) ->.;
    . + guides(color=guide_legend(nrow=2)) -> .;
    . + labs(x="", y="hclust height", color="") ->.;
    ggsave.A4(
        "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.oocyte.hclust.png",
        .,
        width.r=0.45, height.r=0.4
    )
}

## Plot all samples of GSE95477 except for "GSM2514781","GSM2514773"
{
    
    RE.count.info.dt[gse=="GSE95477"] ->.;
    ## remove GSM2514781
    .[SAMPLE %in% c("GSM2514781","GSM2514773")==FALSE] ->.;
    ## bin maternal ages
    .[, maternal.bin.ordered:=factor(c("young", "old")[(as.integer(maternal.age) > 35)+1], levels=c("old", "young")) ]
    ##= prettify cluster
    .[, cluster.renamed.ordered:=factor(c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster])]
    ## prettify stages
    temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    . -> REE.count.boxplot.GSE95477.without.GSM2514781.and.GSM2514773.to.plot.dt
    
    ## start plotting
    REE.count.boxplot.GSE95477.without.GSM2514781.and.GSM2514773.to.plot.dt -> .;
    ggplot(., aes(x=maternal.bin.ordered, y=mean.RE.per.gene, fill=maternal.bin.ordered)) ->.;
    ## add mean comparison
    . + stat_compare_means(comparisons=list(c("old", "young")), method="wilcox.test", method.args=list(alternative="less")) ->.;
    ## plot boxplots
    . + geom_boxplot() ->.;
    ## facet
    . + facet_grid(cluster.renamed.ordered~stage.description.ordered, scales="free_y") ->.;
    ## add theme
    . + theme_pubr() ->.;
    . + theme(axis.text.x = element_blank()) ->.;
    . + scale_y_continuous(limits=c(4.5, 6.5)) ->.;
    . + labs(x="", y="Average # REE-matching edits per gene", fill="Maternal age group") ->.;
    ## save image
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/REE.count.boxplot.GSE95477.without.GSM2514781.and.GSM2514773.png",
        plot=.,
        width.r=0.45, height.r=0.35)
}




{

    copy(REE.count.boxplot.GSE95477.without.GSM2514781.and.GSM2514773.to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='oocyte.GV'] -> A;
    wilcox.test(A[maternal.bin.ordered=="old", mean.RE.per.gene], A[maternal.bin.ordered=="young", mean.RE.per.gene], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="decay at 8-cell", stage="oocyte (GV)", left.group="old", right.group="young", count.of.samples.for.the.left.group=A[maternal.bin.ordered=="old", .N], count.of.samples.for.the.right.group=A[maternal.bin.ordered=="young", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.FILTERED.decay.at.8.cell.oocyte.GV.dt   


    copy(REE.count.boxplot.GSE95477.without.GSM2514781.and.GSM2514773.to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='oocyte.GV'] -> A;
    wilcox.test(A[maternal.bin.ordered=="old", mean.RE.per.gene], A[maternal.bin.ordered=="young", mean.RE.per.gene], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="oocyte (GV)", left.group="old", right.group="young", count.of.samples.for.the.left.group=A[maternal.bin.ordered=="old", .N], count.of.samples.for.the.right.group=A[maternal.bin.ordered=="young", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.FILTERED.others.oocyte.GV.dt   


    copy(REE.count.boxplot.GSE95477.without.GSM2514781.and.GSM2514773.to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='oocyte.MII'] -> A;
    wilcox.test(A[maternal.bin.ordered=="old", mean.RE.per.gene], A[maternal.bin.ordered=="young", mean.RE.per.gene], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="decay at 8-cell", stage="oocyte (MII)", left.group="old", right.group="young", count.of.samples.for.the.left.group=A[maternal.bin.ordered=="old", .N], count.of.samples.for.the.right.group=A[maternal.bin.ordered=="young", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.FILTERED.decay.at.8.cell.oocyte.MII.dt   


    copy(REE.count.boxplot.GSE95477.without.GSM2514781.and.GSM2514773.to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='oocyte.MII'] -> A;
    wilcox.test(A[maternal.bin.ordered=="old", mean.RE.per.gene], A[maternal.bin.ordered=="young", mean.RE.per.gene], alternative="less", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="oocyte (MII)", left.group="old", right.group="young", count.of.samples.for.the.left.group=A[maternal.bin.ordered=="old", .N], count.of.samples.for.the.right.group=A[maternal.bin.ordered=="young", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.FILTERED.others.oocyte.MII.dt   

    list(
        temp.FILTERED.decay.at.8.cell.oocyte.GV.dt,
        temp.FILTERED.others.oocyte.GV.dt,
        temp.FILTERED.decay.at.8.cell.oocyte.MII.dt,
        temp.FILTERED.others.oocyte.MII.dt
    ) -> .;
    rbindlist(., use.names=TRUE) -> .;

    fwrite(., "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/stat.for.REE.count.boxplot.GSE95477.without.GSM2514781.and.GSM2514773.png.csv.gz")
}

