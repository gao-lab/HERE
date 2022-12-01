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

## Plot all 2-cell and 8-cell samples of GSE101571
{
    RE.count.info.dt[gse=="GSE101571"] ->.;
    .[stage %in% c("2-cell", "8-cell")] ->.;
    .[treatment %in% c("amanitin", "control")] ->.;
    ##= prettify cluster
    .[, cluster.renamed.ordered:=factor(c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster])]
    ## prettify stages
    temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## start plotting
    ggplot(., aes(x=treatment, y=mean.RE.per.gene, fill=treatment)) ->.;
    ## add mean comparison
    . + stat_compare_means(comparisons=list(c("amanitin", "control")), method="wilcox.test", method.args=list(alternative="less")) ->.;
    ## plot boxplots
    . + geom_boxplot() ->.;
    ## facet
    . + facet_grid(cluster.renamed.ordered~stage.description.ordered, scales="free_y") ->.;
    ## add theme
    . + theme_pubr() ->.;
    . + theme(axis.text.x = element_blank()) ->.;
    . + scale_y_continuous(expand=expansion(mult=c(0, 0.2)), limits=c(0, 5.5)) ->.;
    . + labs(x="", y="Average # RE-matching edits per gene", fill="Treatment") ->.;
    ## save image
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/RE.count.boxplot.GSE101571.png",
        plot=.,
        width.r=0.45, height.r=0.75)
}
