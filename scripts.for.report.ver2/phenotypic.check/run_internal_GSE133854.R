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

## Plot zygote ~ morula samples of GSE133854
{
    
    RE.count.info.dt[gse=="GSE133854"] ->.;
    .[stage %in% c("zygote", "2-cell", "4-cell", "8-cell", "morula")] ->.;
    ## prettify diseases
    .[disease == "androgenetic", disease.renamed:="AG"]
    .[disease == "parthenogenetic", disease.renamed:="PG"]
    .[is.normal == TRUE & stage != "oocyte.MII", disease.renamed:="BI"]
    .[, disease.renamed.ordered:=factor(disease.renamed, levels=c("AG", "PG", "BI"))]
    ##= prettify cluster
    .[, cluster.renamed.ordered:=factor(c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster])]
    ## prettify stages
    temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    . -> to.plot.dt
    
    ## start plotting
    to.plot.dt -> .;
    ggplot(., aes(x=disease.renamed.ordered, y=mean.RE.per.gene, fill=disease.renamed)) ->.;
    ## add mean comparison
    . + stat_compare_means(comparisons=list(c("AG", "BI"), c("PG", "BI")), method="wilcox.test", method.args=list(alternative="greater")) ->.;
    ## plot boxplots
    . + geom_boxplot() ->.;
    ## facet
    . + facet_grid(cluster.renamed.ordered~stage.description.ordered, scales="free_y") ->.;
    ## add theme
    . + theme_pubr() ->.;
    . + theme(axis.text.x = element_blank()) ->.;
    . + scale_y_continuous(limits=c(0, 9), breaks=seq(0, 8, 2)) ->.;
    . + labs(x="", y="Average # REE-matching edits per gene", fill="Embryo type") ->.;
    ## save image
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/REE.count.boxplot.GSE133854.png",
        plot=.,
        width.r=0.45, height.r=0.75)
}


{

    ## decay at 8-cell; zygote
    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='zygote'] -> A;
    wilcox.test(A[disease.renamed.ordered=="AG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="decay at 8-cell", stage="zygote", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="AG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.decay.at.8.cell.zygote.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='zygote'] -> A;
    wilcox.test(A[disease.renamed.ordered=="PG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="decay at 8-cell", stage="zygote", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="PG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.decay.at.8.cell.zygote.PG.vs.BI.dt   

    ## decay at 8-cell; 2-cell
    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='2-cell'] -> A;
    wilcox.test(A[disease.renamed.ordered=="AG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="decay at 8-cell", stage="2-cell", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="AG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.decay.at.8.cell.2cell.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='2-cell'] -> A;
    wilcox.test(A[disease.renamed.ordered=="PG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="decay at 8-cell", stage="2-cell", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="PG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.decay.at.8.cell.2cell.PG.vs.BI.dt   

    ## decay at 8-cell; 4-cell
    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='4-cell'] -> A;
    wilcox.test(A[disease.renamed.ordered=="AG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="decay at 8-cell", stage="4-cell", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="AG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.decay.at.8.cell.4cell.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='4-cell'] -> A;
    wilcox.test(A[disease.renamed.ordered=="PG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="decay at 8-cell", stage="4-cell", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="PG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.decay.at.8.cell.4cell.PG.vs.BI.dt   


    ## decay at 8-cell; 8-cell
    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='8-cell'] -> A;
    data.table(cluster="decay at 8-cell", stage="8-cell", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="AG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=as.numeric(NA), CI.95pct.lower.bound=as.numeric(NA), CI.95pct.upper.bound=as.numeric(NA)) -> temp.decay.at.8.cell.8cell.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="decay at 8-cell"][stage=='8-cell'] -> A;
    data.table(cluster="decay at 8-cell", stage="8-cell", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="PG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=as.numeric(NA), CI.95pct.lower.bound=as.numeric(NA), CI.95pct.upper.bound=as.numeric(NA)) -> temp.decay.at.8.cell.8cell.PG.vs.BI.dt   

    ## others; zygote
    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='zygote'] -> A;
    wilcox.test(A[disease.renamed.ordered=="AG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="zygote", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="AG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.others.zygote.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='zygote'] -> A;
    wilcox.test(A[disease.renamed.ordered=="PG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="zygote", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="PG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.others.zygote.PG.vs.BI.dt   

    ## others; 2-cell
    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='2-cell'] -> A;
    wilcox.test(A[disease.renamed.ordered=="AG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="2-cell", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="AG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.others.2cell.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='2-cell'] -> A;
    wilcox.test(A[disease.renamed.ordered=="PG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="2-cell", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="PG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.others.2cell.PG.vs.BI.dt   

    ## others; 4-cell
    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='4-cell'] -> A;
    wilcox.test(A[disease.renamed.ordered=="AG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="4-cell", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="AG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.others.4cell.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='4-cell'] -> A;
    wilcox.test(A[disease.renamed.ordered=="PG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="4-cell", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="PG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.others.4cell.PG.vs.BI.dt   


    ## others; 8-cell
    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='8-cell'] -> A;
    wilcox.test(A[disease.renamed.ordered=="AG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="8-cell", left.group="AG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="AG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.others.8cell.AG.vs.BI.dt   

    copy(to.plot.dt) -> A;
    A[cluster.renamed.ordered=="others"][stage=='8-cell'] -> A;
    wilcox.test(A[disease.renamed.ordered=="PG", mean.RE.per.gene], A[disease.renamed.ordered=="BI", mean.RE.per.gene], alternative="greater", paired=FALSE, conf.int=TRUE) -> A.wilcox
    data.table(cluster="others", stage="8-cell", left.group="PG", right.group="BI", count.of.samples.for.the.left.group=A[disease.renamed.ordered=="PG", .N], count.of.samples.for.the.right.group=A[disease.renamed.ordered=="BI", .N], difference.in.location=A.wilcox$estimate, CI.95pct.lower.bound=A.wilcox$conf.int[1], CI.95pct.upper.bound=A.wilcox$conf.int[2]) -> temp.others.8cell.PG.vs.BI.dt   


    list(
        temp.decay.at.8.cell.zygote.AG.vs.BI.dt,   
        temp.decay.at.8.cell.zygote.PG.vs.BI.dt,   
        temp.decay.at.8.cell.2cell.AG.vs.BI.dt,   
        temp.decay.at.8.cell.2cell.PG.vs.BI.dt,   
        temp.decay.at.8.cell.4cell.AG.vs.BI.dt,   
        temp.decay.at.8.cell.4cell.PG.vs.BI.dt,
        temp.decay.at.8.cell.8cell.AG.vs.BI.dt,   
        temp.decay.at.8.cell.8cell.PG.vs.BI.dt,
        temp.others.zygote.AG.vs.BI.dt,   
        temp.others.zygote.PG.vs.BI.dt,   
        temp.others.2cell.AG.vs.BI.dt,   
        temp.others.2cell.PG.vs.BI.dt,   
        temp.others.4cell.AG.vs.BI.dt,   
        temp.others.4cell.PG.vs.BI.dt,
        temp.others.8cell.AG.vs.BI.dt,   
        temp.others.8cell.PG.vs.BI.dt
    ) -> .;
    rbindlist(., use.names=TRUE) -> .;

    fwrite(., "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/stat.for.REE.count.boxplot.GSE133854.png.csv.gz")
}
