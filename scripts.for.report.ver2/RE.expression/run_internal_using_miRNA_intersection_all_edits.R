library("data.table")
library("readxl")
library("ggalluvial")
library("foreach")
library("ggtext")
library("glue")
library("statpsych")
source("./scripts/common/ggpubr.A4.R")

## We focus on normal early stages only.

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.and.expression.using.miRNA.intersection.all.edits"
dir.create(output.directory, recursive=TRUE)

{
    c(
        "oocyte.GV->oocyte.MII"="oocyte (GV)\n-> oocyte (MII)",
        "oocyte.MII->zygote"="oocyte (MII)\n-> zygote",
        "zygote->2-cell"="zygote\n-> 2-cell",
        "2-cell->4-cell"="2-cell\n-> 4-cell",
        "4-cell->8-cell"="4-cell\n-> 8-cell"
    ) -> .;
} -> stage.pair.named.vector



## 1. Get % edited transcripts per gene per sample
{
    ##
    ## 1.1. get editing level for each edited events on valid genes in normal early stages
    {
        fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz") -> .;
        ## keep edited events only
        .[is.na(AC) == FALSE] -> .;
        ## keep normal early stages only
        .[is.normal == TRUE & stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell')] -> .;
        ## .[, .N, list(SAMPLE, ID, Gene_ID, Annotation)][N>1]
        ## keep relevant columns and take the unique value
        ## (reason of duplicated row: each single Annotation is written separately for a given edit that has multiple annotations on a single gene)
        ## .[, .N, list(SAMPLE, ID, Gene_ID, Annotation)][N>1] : has 0 rows, but will have many rows if Annotation is removed
        unique(.[, list(
            SAMPLE, gse, stage,
            Gene_ID, Gene_Name,
            ID, CHROM, POS, AC, AN, AF
        )]) -> .;
    } -> ..RE.matching.edits.edited.only.valid.genes.only.in.normal.early.stages.only.dt;
    ##
} -> RE.matching.edits.edited.only.valid.genes.only.in.normal.early.stages.only.dt;

## 2. Get FPKM per gene per sample
{
    ##
    ## 2.1. get sample info of all normal early stage samples
    {
        fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt") -> .;
        ## keep normal early stage samples
        .[is.normal==TRUE & stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell')] -> .;
    } -> ..normal.early.stage.samples.info.dt
    ##
    ## 2.2. get FPKM for all normal early stage samples (only samples in the 2,071 set are considered; the rest are discarded due to insufficient read length (e.g, GSM1160130, is 8-cell, but has read length of 49*2 < 75*2)
    {
        read.table("result/BS06_1__get_expression_level/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/combined.gexpr.FPKM.matrix.txt") -> .;
        ## keep matched samples only
        .[, intersect(colnames(.), ..normal.early.stage.samples.info.dt[, gsm])] -> .;
        ## transform to 3-way table
        data.table(Gene_ID=rownames(.), .) -> .;
        melt(., id.vars="Gene_ID", variable.name="SAMPLE", value.name="FPKM") -> .;
    } -> ..combined.gexpr.FPKM.melt.normal.early.stage.samples.only.dt;
    ##
    ..combined.gexpr.FPKM.melt.normal.early.stage.samples.only.dt
} -> combined.gexpr.FPKM.melt.normal.early.stage.samples.only.dt

## 3. Merge editing level and FPKM per (edit x gene), and compute the spearman correlation between them (within sample)
{
    ## 3.1. merge the two tables (within-sample style)
    {
        merge(x=RE.matching.edits.edited.only.valid.genes.only.in.normal.early.stages.only.dt,
              y=combined.gexpr.FPKM.melt.normal.early.stage.samples.only.dt,
              by.x=c("SAMPLE", "Gene_ID"), by.y=c("SAMPLE", "Gene_ID"),
              all.x=TRUE, all.y=FALSE) -> .;
    } -> ..RE.matching.edits.edited.only.and.FPKM.valid.genes.only.in.normal.early.stages.only.dt;
    ##
    ## 3.2. compute the spearman correlation (within samples)
    {
        ..RE.matching.edits.edited.only.and.FPKM.valid.genes.only.in.normal.early.stages.only.dt -> .;
        .[, list(
            spearman.cor.within.sample=cor(AF, FPKM, method="spearman")
        ), list(SAMPLE, gse, stage)] -> .;
    } -> ..spearman.cor.within.sample.across.genes.dt
    fwrite(..spearman.cor.within.sample.across.genes.dt, glue("{output.directory}/spearman.cor.within.sample.across.genes.dt.gz"))
} -> spearman.cor.within.sample.across.genes.dt


## 4. Merge editing level and FPKM per (edit x gene), and compute the spearman correlation between them (consecutive-stage style)
{
    ..temp.list <- list()
    ##
    ## 4.1. construct consecutive-stage table
    {
        data.table(
            stage.previous=c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell"),
            stage.current=c("oocyte.MII", "zygote", "2-cell", "4-cell", "8-cell")
        ) -> .;
    } -> ..consecutive.stage.dt
    ..temp.list[["consecutive.stage.dt"]] <- ..consecutive.stage.dt
    ##
    ## 4.2. compute median AF per (gene, stage)
    {
        RE.matching.edits.edited.only.valid.genes.only.in.normal.early.stages.only.dt -> .;
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
        merge(x=combined.gexpr.FPKM.melt.normal.early.stage.samples.only.dt,
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
        merge(x=..consecutive.stage.dt, y=..median.AF.per.gene.and.stage.dt[, list(stage, Gene_ID, Gene_Name, ID, CHROM, POS, AF.median.in.stage.previous=AF.median)],
              by.x="stage.previous", by.y="stage",
              all.x=TRUE, all.y=FALSE) -> .;
        merge(x=., y=..median.FPKM.per.gene.and.stage.dt[, list(stage, Gene_ID, FPKM.median.in.stage.current=FPKM.median)],
              by.x=c("stage.current", "Gene_ID"), by.y=c("stage", "Gene_ID"),
              all=TRUE, all.y=FALSE) -> .;
    } -> ..median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt
    fwrite(..median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt, glue("{output.directory}/median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt.gz"))
    ..temp.list[["median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt"]] <- copy(..median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt)
    ##
    ..temp.list
} -> temp.list
consecutive.stage.dt <- temp.list[["consecutive.stage.dt"]]
median.AF.per.gene.and.stage.dt <- temp.list[["median.AF.per.gene.and.stage.dt"]]
median.FPKM.per.gene.and.stage.dt <- temp.list[["median.FPKM.per.gene.and.stage.dt"]]
median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt <- temp.list[["median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt"]]


## 5. Prepare for annotations with maternal info and MBS info
{
    ..temp.list <- list()
    
    ## 5.1. Get maternal gene info
    {
        fread(
            "result/S42_1__annotate_embryonic_genes/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/201221-fifth-phenotype-collection/combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt.txt.gz"
        ) -> .;
        .[
            cluster %in% c("maternal.decay", "maternal.others"),
            list(Gene_ID=gene.id, cluster)
        ] -> .;
    } -> ..maternal.gene.info.dt;
    ..temp.list[["maternal.gene.info.dt"]] <- ..maternal.gene.info.dt
    
    ## 5.2. Annotate MBS type per edit
    {
        fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA.intersection.all.edits/all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.summary.gene.and.edit.level.dt.gz") -> .;
    } -> ..MBS.annotation.per.edit.and.gene.dt;
    ..temp.list[["MBS.annotation.per.edit.and.gene.dt"]] <- ..MBS.annotation.per.edit.and.gene.dt
    
    ..temp.list
} -> temp.list
maternal.gene.info.dt <- temp.list[["maternal.gene.info.dt"]]
MBS.annotation.per.edit.and.gene.dt <- temp.list[["MBS.annotation.per.edit.and.gene.dt"]]

## 6. Annotate maternal gene info and MBS annotation onto median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt
{
    merge(x=median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.dt, y=maternal.gene.info.dt,
          by.x="Gene_ID", by.y="Gene_ID",
          all.x=TRUE, all.y=FALSE) ->.;
    merge(x=., y=MBS.annotation.per.edit.and.gene.dt,
          by.x=c("Gene_ID", "POS"), by.y=c("gene.id", "edit.POS"),
          all.x=TRUE, all.y=FALSE) -> .;
    . -> ..median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.with.maternal.cluster.and.MBS.annotation.dt
    fwrite(..median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.with.maternal.cluster.and.MBS.annotation.dt, glue("{output.directory}/median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.with.maternal.cluster.and.MBS.annotation.dt.gz"))
    ..median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.with.maternal.cluster.and.MBS.annotation.dt
} -> median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.with.maternal.cluster.and.MBS.annotation.dt

## 7. Plot scatterplots for AF.median ~ FPKM.median under different conditions
{
    
    ## 7.1. plot all gene x edit pairs
    {
        
        copy(median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.with.maternal.cluster.and.MBS.annotation.dt) -> .;
        .[, list(
            spearman.n=.N, 
            spearman.cor=cor(AF.median.in.stage.previous, FPKM.median.in.stage.current, method="spearman"),
            spearman.cor.test.pvalue=cor.test(AF.median.in.stage.previous, FPKM.median.in.stage.current, method="spearman", alternative="two.sided")$p.value,
            spearman.cor.test.95CI.LL=ci.spear(alpha=0.5, y=AF.median.in.stage.previous, x=FPKM.median.in.stage.current)[, "LL"],
            spearman.cor.test.95CI.UL=ci.spear(alpha=0.5, y=AF.median.in.stage.previous, x=FPKM.median.in.stage.current)[, "UL"]
        ), list(stage.previous, stage.current)] -> ..temp.cor.dt;
        ..temp.cor.dt[, spearman.cor.test.pvalue.adjusted:=p.adjust(spearman.cor.test.pvalue, method="BH")] -> ..temp.cor.dt;
        ..temp.cor.dt[spearman.cor.test.pvalue.adjusted >= 2.2e-16, spearman.cor.test.pvalue.adjusted.prettified:=format(spearman.cor.test.pvalue.adjusted, scientific=TRUE)] -> ..temp.cor.dt;
        ..temp.cor.dt[spearman.cor.test.pvalue.adjusted < 2.2e-16, spearman.cor.test.pvalue.adjusted.prettified:="<2.2e-16"] -> ..temp.cor.dt;
        merge(x=., y=..temp.cor.dt, by=c("stage.previous", "stage.current"), all.x=TRUE, all.y=FALSE) -> .;
        .[, stage.pair.ordered:=factor(stage.pair.named.vector[paste(sep="", stage.previous, "->", stage.current)], levels=stage.pair.named.vector)] -> ..to.plot.dt;
        
        ggplot(..to.plot.dt[FPKM.median.in.stage.current > 1e-3], aes(x=AF.median.in.stage.previous, y=FPKM.median.in.stage.current)) -> .;
        . + geom_point(size = 1, alpha=0.3, shape=4) -> .;
        ##. + geom_density_2d_filled(alpha=0.5) -> .;
        . + geom_smooth(method=lm, color="#EB9D35") -> .;
        . + geom_text(data=unique(..to.plot.dt[, list(stage.pair.ordered, spearman.cor.test.pvalue.adjusted.prettified, spearman.cor.test.pvalue.adjusted, spearman.n, spearman.cor, spearman.cor.test.95CI.LL, spearman.cor.test.95CI.UL, AF.median.in.stage.previous=0.5, FPKM.median.in.stage.current=2e-3)]), mapping=aes(label=glue("Count of genes used:\n{spearman.n}\nSpearman's rho:\n{round(spearman.cor, 4)}\nBH-adjusted p-value:\n{spearman.cor.test.pvalue.adjusted.prettified}\n95% CI:\n[{round(spearman.cor.test.95CI.LL, 4)}, {round(spearman.cor.test.95CI.UL, 4)}]")), size=3) -> .;
        . + scale_x_continuous(limits=c(0.1, 1), breaks=c(0.1, 0.4, 0.7, 1)) -> .;
        . + scale_y_log10(limits=c(1e-3, 1e4), breaks=c(1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4), labels=function(x) {sub(pattern="(\\.0+|0+)$", replacement="", x=format(x, scientific=FALSE))}) -> .;
        . + theme_pubr(base_size=10) -> .;
        . + theme(axis.text.x=element_text(angle=45, hjust=1)) -> .;
        . + facet_grid(~stage.pair.ordered) -> .;
        . + labs(x="median(editing level) in the previous stage\n(>=0.1 due to identification pipeline)", y="median(FPKM) in the current stage") -> .;
        ggsave.A4(filename=glue("{output.directory}/AF.median.vs.FPKM.median.scatterplot.overall.png"), plot=., width.r=1, height.r=0.8)

        
    }
    
    ## 7.2. plot all gene x edit pairs, targets of maternal clearance or not
    {
        copy(median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.with.maternal.cluster.and.MBS.annotation.dt) -> .;
        .[cluster %in% c("maternal.decay", "maternal.others")] -> .;
        .[, spearman.cor:=cor(AF.median.in.stage.previous, FPKM.median.in.stage.current, method="spearman"), list(stage.previous, stage.current, cluster)] -> .;
        .[, stage.pair.ordered:=factor(stage.pair.named.vector[paste(sep="", stage.previous, "->", stage.current)], levels=stage.pair.named.vector)] -> ..to.plot.dt;
        ggplot(..to.plot.dt[FPKM.median.in.stage.current > 1e-3], aes(x=AF.median.in.stage.previous, y=FPKM.median.in.stage.current)) -> .;
        . + geom_point(size = 1, alpha=0.3, shape=4) -> .;
        ##. + geom_density_2d_filled(alpha=0.5) -> .;
        . + geom_smooth(method=lm, color="#EB9D35") -> .;
        . + geom_text(data=unique(..to.plot.dt[, list(stage.pair.ordered, cluster, spearman.cor, AF.median.in.stage.previous=0.5, FPKM.median.in.stage.current=2e-3)]), mapping=aes(label=glue("Spearman's rho:\n{round(spearman.cor, 4)}")), size=3) -> .;
        . + scale_x_continuous(limits=c(0.1, 1), breaks=c(0.1, 0.4, 0.7, 1)) -> .;
        . + scale_y_log10(limits=c(1e-3, 1e4), breaks=c(1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4), labels=function(x) {sub(pattern="(\\.0+|0+)$", replacement="", x=format(x, scientific=FALSE))}) -> .;
        . + theme_pubr(base_size=10) -> .;
        . + theme(axis.text.x=element_text(angle=45, hjust=1)) -> .;
        . + facet_grid(cluster~stage.pair.ordered) -> .;
        . + labs(x="median(editing level) in the previous stage\n(>=0.1 due to identification pipeline)", y="median(FPKM) in the current stage") -> .;
        ggsave.A4(filename=glue("{output.directory}/AF.median.vs.FPKM.median.scatterplot.wrt.cluster.png"), plot=., width.r=0.9, height.r=0.4)
    }

    ## 7.3. plot all gene x edit pairs, targets of maternal clearance or not, MBS annotation type
    {
        copy(median.AF.stage.previous.and.median.FPKM.stage.current.per.gene.with.maternal.cluster.and.MBS.annotation.dt) -> .;
        .[cluster %in% c("maternal.decay", "maternal.others")] -> .;
        .[, spearman.cor:=cor(AF.median.in.stage.previous, FPKM.median.in.stage.current, method="spearman"), list(stage.previous, stage.current, cluster, gains.miRNA.sites)] -> .;
        .[, stage.pair.ordered:=factor(stage.pair.named.vector[paste(sep="", stage.previous, "->", stage.current)], levels=stage.pair.named.vector)] -> .;
        .[, gains.miRNA.sites.prettified:=c("REE is not\nMBS-gaining", "REE is\nMBS-gaining")[(gains.miRNA.sites+1)]][is.na(gains.miRNA.sites.prettified) == TRUE, gains.miRNA.sites.prettified:="no MBS identified\nregardless of REE editing"] -> .;
        .[, cluster.prettified:=c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster]] -> .;
        . -> ..to.plot.dt;
        ggplot(..to.plot.dt[FPKM.median.in.stage.current > 1e-3], aes(x=AF.median.in.stage.previous, y=FPKM.median.in.stage.current)) -> .;
        . + geom_point(size = 1, alpha=0.3, shape=4) -> .;
        ##. + geom_density_2d_filled(alpha=0.5) -> .;
        . + geom_smooth(method=lm, color="#EB9D35") -> .;
        . + geom_text(data=unique(..to.plot.dt[, list(stage.pair.ordered, cluster.prettified, gains.miRNA.sites.prettified, spearman.cor, AF.median.in.stage.previous=0.5, FPKM.median.in.stage.current=2e-3)]), mapping=aes(label=glue("Spearman's rho:\n{format(round(spearman.cor, 4), scientific=FALSE)}")), size=3) -> .;
        . + scale_x_continuous(limits=c(0.1, 1), breaks=c(0.1, 0.4, 0.7, 1)) -> .;
        . + scale_y_log10(limits=c(1e-3, 1e4), breaks=c(1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4), labels=function(x) {sub(pattern="(\\.0+|0+)$", replacement="", x=format(x, scientific=FALSE))}) -> .;
        . + theme_pubr(base_size=10) -> .;
        . + theme(axis.text.x=element_text(angle=45, hjust=1)) -> .;
        . + facet_grid(cluster.prettified + gains.miRNA.sites.prettified ~stage.pair.ordered) -> .;
        . + labs(x="median(editing level) in the previous stage\n(>=0.1 due to identification pipeline)", y="median(FPKM) in the current stage") -> .;
        ggsave.A4(filename=glue("{output.directory}/REE.AF.median.vs.FPKM.median.scatterplot.wrt.cluster.and.MBS.annotation.png"), plot=., width.r=0.9, height.r=1)
    }
    
}

## [SAMPLE, stage, Gene, RE edit (ID, editing level)] 
## |
## SAMPLE, stage, Gene, % RE-edited transcripts of this gene (lower bound-upper bound)
## | <- [SAMPLE, Gene, FPKM]
## SAMPLE, stage, Gene, % RE-edited transcripts of this gene (lower bound-upper bound), FPKM
## |
## stage-before, stage-after, spearman(% RE-edited transcripts, FPKM | across all genes (could be classified into maternal clearance targets and others)) (should be overall negative if RE helps to decay)




## 8. Merge editing level diff and FPKM diff per (edit x gene), and compute the spearman correlation between them (can only be consecutive-stage style)
{
    ..temp.list <- list()
    
    ## 8.1. compute diff of median AF per (gene, stage)
    {
        ## merge previous stage RE
        merge(x=consecutive.stage.dt,
              y=median.AF.per.gene.and.stage.dt[, list(stage, Gene_ID, Gene_Name, ID, CHROM, POS, AF.median.in.stage.previous=AF.median)],
              by.x="stage.previous", by.y="stage",
              all.x=TRUE, all.y=FALSE) -> .;
        ## merge current stage RE
        merge(x=.,
              y=median.AF.per.gene.and.stage.dt[, list(stage, Gene_ID, Gene_Name, ID, CHROM, POS, AF.median.in.stage.current=AF.median)],
              by.x=c("stage.current", "Gene_ID", "Gene_Name", "ID", "CHROM", "POS"),
              by.y=c("stage", "Gene_ID", "Gene_Name", "ID", "CHROM", "POS"),
              all.x=TRUE, all.y=FALSE) -> .;
        ## label two types of transitions (RE.to.RE and RE.to.non-RE)
        .[is.na(AF.median.in.stage.current) == FALSE, transition.type:="RE.to.RE"] -> .;
        .[is.na(AF.median.in.stage.current) == TRUE, transition.type:="RE.to.non-RE"] -> .;
        ## compute AF diff for RE.to.RE
        .[transition.type == "RE.to.RE", diff.AF.median := AF.median.in.stage.current - AF.median.in.stage.previous] -> .;
    } -> ..diff.median.AF.per.gene.and.stage.pair.dt
    ..temp.list[["diff.median.AF.per.gene.and.stage.pair.dt"]] <- ..diff.median.AF.per.gene.and.stage.pair.dt
   
    ## 8.2. compute diff of median FPKM per (gene, stage pair)
    {
        ## merge previous stage FPKM
        merge(x=consecutive.stage.dt,
              y=median.FPKM.per.gene.and.stage.dt[, list(stage, Gene_ID, FPKM.median.in.stage.previous=FPKM.median)],
              by.x="stage.previous", by.y="stage",
              all.x=TRUE, all.y=FALSE) -> .;
        ## merge current stage FPKM
        merge(x=.,
              y=median.FPKM.per.gene.and.stage.dt[, list(stage, Gene_ID, FPKM.median.in.stage.current=FPKM.median)],
              by.x=c("stage.current", "Gene_ID"),
              by.y=c("stage", "Gene_ID"),
              all.x=TRUE, all.y=FALSE) -> .;        
        ## compute FPKM diff
        .[, diff.FPKM.median := FPKM.median.in.stage.current - FPKM.median.in.stage.previous] -> .;
        .[, diff.log2.corrected.FPKM.median := log2((FPKM.median.in.stage.current + 1e-4)/(FPKM.median.in.stage.previous + 1e-4))] -> .;
    } -> ..diff.median.FPKM.per.gene.and.stage.pair.dt
    ..temp.list[["diff.median.FPKM.per.gene.and.stage.pair.dt"]] <- ..diff.median.FPKM.per.gene.and.stage.pair.dt

    ## 8.3. merge the two tables
    {
        merge(x=..diff.median.AF.per.gene.and.stage.pair.dt,
              y=..diff.median.FPKM.per.gene.and.stage.pair.dt,
              by.x=c("stage.previous", "stage.current", "Gene_ID"),
              by.y=c("stage.previous", "stage.current", "Gene_ID"),
              all.x=TRUE, all.y=FALSE) -> .;
    } -> ..diff.median.AF.and.FPKM.per.gene.and.stage.pair.dt
    fwrite(..diff.median.AF.and.FPKM.per.gene.and.stage.pair.dt, glue("{output.directory}/diff.median.AF.and.FPKM.per.gene.and.stage.pair.dt.gz"))
    ..temp.list[["diff.median.AF.and.FPKM.per.gene.and.stage.pair.dt"]] <- copy(..diff.median.AF.and.FPKM.per.gene.and.stage.pair.dt)
    ..temp.list
} -> temp.list
diff.median.AF.and.FPKM.per.gene.and.stage.pair.dt <- temp.list[["diff.median.AF.and.FPKM.per.gene.and.stage.pair.dt"]]

## 9. Annotate maternal gene info and MBS annotation onto diff.median.AF.and.FPKM.per.gene.and.stage.pair.dt
{
    merge(x=diff.median.AF.and.FPKM.per.gene.and.stage.pair.dt, y=maternal.gene.info.dt,
          by.x="Gene_ID", by.y="Gene_ID",
          all.x=TRUE, all.y=FALSE) ->.;
    merge(x=., y=MBS.annotation.per.edit.and.gene.dt,
          by.x=c("Gene_ID", "POS"), by.y=c("gene.id", "edit.POS"),
          all.x=TRUE, all.y=FALSE) -> .;
    . -> ..diff.median.AF.and.FPKM.per.gene.and.stage.pair.with.maternal.cluster.and.MBS.annotation.dt
    fwrite(..diff.median.AF.and.FPKM.per.gene.and.stage.pair.with.maternal.cluster.and.MBS.annotation.dt, glue("{output.directory}/diff.median.AF.and.FPKM.per.gene.and.stage.pair.with.maternal.cluster.and.MBS.annotation.dt.gz"))
    ..diff.median.AF.and.FPKM.per.gene.and.stage.pair.with.maternal.cluster.and.MBS.annotation.dt
} -> diff.median.AF.and.FPKM.per.gene.and.stage.pair.with.maternal.cluster.and.MBS.annotation.dt


## 10. Plot scatterplots for diff.AF.median ~ diff.FPKM.median under different conditions
{
        
    ## 10.1. plot histogram of transition types, focusing on targets of maternal clearance only, focusing on cleared genes only
    {
        ##
        copy(diff.median.AF.and.FPKM.per.gene.and.stage.pair.with.maternal.cluster.and.MBS.annotation.dt) -> .;
        .[cluster %in% c("maternal.decay", "maternal.others")] -> .;
        .[diff.FPKM.median < 0] -> .;
        .[, stage.pair.ordered:=factor(stage.pair.named.vector[paste(sep="", stage.previous, "->", stage.current)], levels=stage.pair.named.vector)] -> .;
        .[, transition.type.prettified:=c("RE.to.RE"="is still a REE in the latter stage", "RE.to.non-RE"="is not a REE in the latter stage")[transition.type]] -> .;
        .[, cluster.prettified:=c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster]] -> .;
        . -> ..to.plot.dt;
        ##
        ggplot(..to.plot.dt, aes(x=transition.type.prettified, fill=transition.type.prettified)) -> .;
        . + stat_count(position="dodge") -> .;
        . + theme_pubr(base_size=10) -> .;
        . + theme(axis.text.x=element_blank()) -> .;
        . + facet_grid(cluster.prettified~stage.pair.ordered) -> .;
        . + labs(x="", y="# (gene, REE) pairs\nwhose median gene FPKM drops\nupon stage transition", fill="REE group") -> .;
        ggsave.A4(filename=glue("{output.directory}/REE.transition.type.histogram.gene.median.FPKM.drop.only.wrt.cluster.png"), plot=., width.r=0.75, height.r=0.3)
        ##
    } -> not.used.variable
        

    ## 10.2. plot all gene x edit pairs, focusing on targets of maternal clearance only, focusing on cleared genes only
    {
        
        copy(diff.median.AF.and.FPKM.per.gene.and.stage.pair.with.maternal.cluster.and.MBS.annotation.dt) -> .;
        .[cluster %in% c("maternal.decay", "maternal.others")] -> .;
        .[diff.FPKM.median < 0] -> .;
        .[transition.type == 'RE.to.RE'] -> .;
        .[, stage.pair.ordered:=factor(stage.pair.named.vector[paste(sep="", stage.previous, "->", stage.current)], levels=stage.pair.named.vector)] -> .;
        .[, cluster.prettified:=c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster]] -> .;
        . -> ..to.plot.dt;
        
        ggplot(..to.plot.dt, aes(x=cluster.prettified, y=diff.AF.median, fill=cluster.prettified)) -> .;
        . + geom_boxplot() -> .;
        . + geom_hline(yintercept=0, linetype="dashed") -> .;
        . + stat_compare_means(comparisons=list(c("decay at 8-cell", "others")), method="wilcox.test", method.args=list(alternative="less")) ->.;
        . + theme_pubr(base_size=10) -> .;
        . + theme(axis.text.x=element_blank()) -> .;
        . + facet_grid(~stage.pair.ordered, scales="free_y") -> .;
        . + labs(x="", y="median(editing level) in the latter stage\nminus median(editing level) in the former stage\non genes with median FPKM dropped upon stage transition", fill="") -> .;
        ggsave.A4(filename=glue("{output.directory}/diff.AF.median.boxplot.gene.median.FPKM.drop.only.wrt.cluster.png"), plot=., width.r=0.7, height.r=0.5)
        
    } -> not.used.variable


}
