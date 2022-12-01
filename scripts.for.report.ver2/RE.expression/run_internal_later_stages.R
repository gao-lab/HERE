library("data.table")
library("readxl")
library("ggalluvial")
library("foreach")
library("ggtext")
library("glue")
source("./scripts/common/ggpubr.A4.R")

## We focus on normal cross stages only.

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.and.expression.later.stages"
dir.create(output.directory, recursive=TRUE)

{
    c(
        "8-cell->morula"="8-cell\n-> morula",
        "8-cell->blastocyst.late"="8-cell\n-> blastocyst (late)",
        "8-cell->ICM"="8-cell\n-> ICM"
    ) -> .;
} -> stage.pair.named.vector



## 1. Get % edited transcripts per gene per sample
{
    ##
    ## 1.1. get editing level for each edited events on valid genes in normal cross stages
    {
        fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz") -> .;
        ## keep edited events only
        .[is.na(AC) == FALSE] -> .;
        ## keep normal cross stages only
        .[is.normal == TRUE & stage %in% c("8-cell", "morula", "blastocyst.late", "ICM")] -> .;
        ## .[, .N, list(SAMPLE, ID, Gene_ID, Annotation)][N>1]
        ## keep relevant columns and take the unique value
        ## (reason of duplicated row: each single Annotation is written separately for a given edit that has multiple annotations on a single gene)
        ## .[, .N, list(SAMPLE, ID, Gene_ID, Annotation)][N>1] : has 0 rows, but will have many rows if Annotation is removed
        unique(.[, list(
            SAMPLE, gse, stage,
            Gene_ID, Gene_Name,
            ID, CHROM, POS, AC, AN, AF
        )]) -> .;
    } -> ..RE.matching.edits.edited.only.valid.genes.only.in.normal.cross.stages.only.dt;
    ##
} -> RE.matching.edits.edited.only.valid.genes.only.in.normal.cross.stages.only.dt;


## 2. Get FPKM per gene per sample
{
    ##
    ## 2.1. get sample info of all normal cross stage samples
    {
        fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt") -> .;
        ## keep normal cross stage samples
        .[is.normal==TRUE & stage %in% c("8-cell", "morula", "blastocyst.late", "ICM")] -> .;
    } -> ..normal.cross.stage.samples.info.dt    
    ##
    ## 2.2. get FPKM for all normal cross stage samples (only samples in the 2,071 set are considered; the rest are discarded due to insufficient read length (e.g, GSM1160130, is 8-cell, but has read length of 49*2 < 75*2)
    {
        read.table("result/BS06_1__get_expression_level/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/combined.gexpr.FPKM.matrix.txt") -> .;
        ## keep matched samples only
        .[, intersect(colnames(.), ..normal.cross.stage.samples.info.dt[, gsm])] -> .;
        ## transform to 3-way table
        data.table(Gene_ID=rownames(.), .) -> .;
        melt(., id.vars="Gene_ID", variable.name="SAMPLE", value.name="FPKM") -> .;
    } -> ..combined.gexpr.FPKM.melt.normal.cross.stage.samples.only.dt;
    ##
    fwrite(..combined.gexpr.FPKM.melt.normal.cross.stage.samples.only.dt, glue("{output.directory}/combined.gexpr.FPKM.melt.normal.cross.stage.samples.only.dt.gz"))
    ..combined.gexpr.FPKM.melt.normal.cross.stage.samples.only.dt
} -> combined.gexpr.FPKM.melt.normal.cross.stage.samples.only.dt




## 4. Merge editing level and FPKM per (edit x gene), and compute the spearman correlation between them (cross-stage style)
{
    
    ..temp.list <- list()
    ##
    ## 4.1. construct cross-stage table
    {
        data.table(
            stage.previous=c("8-cell", "8-cell", "8-cell"),
            stage.current=c("morula", "blastocyst.late", "ICM")
        ) -> .;
    } -> ..cross.stage.dt
    ..temp.list[["cross.stage.dt"]] <- ..cross.stage.dt
    
    ##
    ## 4.2. compute median AF per (gene, stage)
    {
        RE.matching.edits.edited.only.valid.genes.only.in.normal.cross.stages.only.dt -> .;
        .[, list(AF.median=median(AF)),
          list(
              stage,
              Gene_ID, Gene_Name,
              ID, CHROM, POS)] -> .;
    } -> ..median.AF.per.gene.and.stage.dt
    fwrite(..median.AF.per.gene.and.stage.dt, glue("{output.directory}/median.AF.per.gene.and.stage.dt.gz"))
    ..temp.list[["median.AF.per.gene.and.stage.dt"]] <- ..median.AF.per.gene.and.stage.dt
    
    ##    
    ## 4.3. compute median FPKM per (gene, stage)
    {
        ## read stage info first
        fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt") -> .;
        ## append stage info to FPKM table
        merge(x=combined.gexpr.FPKM.melt.normal.cross.stage.samples.only.dt,
              y=.[, list(gsm, stage)],
              by.x="SAMPLE", by.y="gsm",
              all.x=TRUE, all.y=FALSE
              )-> .;
        ## compute FPKM median
        .[, list(FPKM.median=median(FPKM)), list(stage, Gene_ID)] -> ..temp.FPKM.dt;
    } -> ..median.FPKM.per.gene.and.stage.dt;
    fwrite(..median.FPKM.per.gene.and.stage.dt, glue("{output.directory}/median.FPKM.per.gene.and.stage.dt.gz"))
    ..temp.list[["median.FPKM.per.gene.and.stage.dt"]] <- ..median.FPKM.per.gene.and.stage.dt
    
    ##
    ## 4.4. merge the three tables
    {
        merge(x=..cross.stage.dt, y=..median.AF.per.gene.and.stage.dt[, list(stage, Gene_ID, Gene_Name, ID, CHROM, POS, AF.median.in.stage.previous=AF.median)],
              by.x="stage.previous", by.y="stage",
              all.x=TRUE, all.y=FALSE, allow.cartesian=TRUE) -> .; ## use allow.cartesian=TRUE because the 'stage.previous' has duplicated 8-cells (which is what we'd like to examine)
        merge(x=., y=..median.FPKM.per.gene.and.stage.dt[, list(stage, Gene_ID, FPKM.median.in.stage.current=FPKM.median)],
              by.x=c("stage.current", "Gene_ID"), by.y=c("stage", "Gene_ID"),
              all=TRUE, all.y=FALSE) -> .;
    } -> ..median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt
    fwrite(..median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt, glue("{output.directory}/median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt.gz"))
    ..temp.list[["median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt"]] <- copy(..median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt)
    
    ##
    ..temp.list
} -> temp.list
cross.stage.dt <- temp.list[["cross.stage.dt"]]
median.AF.per.gene.and.stage.dt <- temp.list[["median.AF.per.gene.and.stage.dt"]]
median.FPKM.per.gene.and.stage.dt <- temp.list[["median.FPKM.per.gene.and.stage.dt"]]
median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt <- temp.list[["median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt"]]



## 7. Plot scatterplots for AF.median ~ FPKM.median under different conditions
{
    
    ## 7.1. plot all gene x edit pairs
    {
        copy(median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt) -> .;
        .[, spearman.cor:=cor(AF.median.in.stage.previous, FPKM.median.in.stage.current, method="spearman"), list(stage.previous, stage.current)] -> .;
        .[, stage.pair.ordered:=factor(stage.pair.named.vector[paste(sep="", stage.previous, "->", stage.current)], levels=stage.pair.named.vector)] -> ..to.plot.dt;
        ggplot(..to.plot.dt[FPKM.median.in.stage.current > 1e-3], aes(x=AF.median.in.stage.previous, y=FPKM.median.in.stage.current)) -> .;
        . + geom_point(size = 1, alpha=0.3, shape=4) -> .;
        ##. + geom_density_2d_filled(alpha=0.5) -> .;
        . + geom_smooth(method=lm, color="#EB9D35") -> .;
        . + geom_text(data=unique(..to.plot.dt[, list(stage.pair.ordered, spearman.cor, AF.median.in.stage.previous=0.5, FPKM.median.in.stage.current=2e-3)]), mapping=aes(label=glue("Spearman's rho:\n{round(spearman.cor, 4)}")), size=3) -> .;
        . + scale_x_continuous(limits=c(0.1, 1), breaks=c(0.1, 0.4, 0.7, 1)) -> .;
        . + scale_y_log10(limits=c(1e-3, 1e4), breaks=c(1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4), labels=function(x) {sub(pattern="(\\.0+|0+)$", replacement="", x=format(x, scientific=FALSE))}) -> .;
        . + theme_pubr(base_size=10) -> .;
        . + theme(axis.text.x=element_text(angle=45, hjust=1)) -> .;
        . + facet_grid(~stage.pair.ordered) -> .;
        . + labs(x="median(editing level) in the previous stage\n(>=0.1 due to identification pipeline)", y="median(FPKM) in the current stage") -> .;
        ggsave.A4(filename=glue("{output.directory}/AF.median.vs.FPKM.median.scatterplot.overall.png"), plot=., width.r=0.9, height.r=0.3)
    }
    
}

