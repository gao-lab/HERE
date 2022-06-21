library("data.table")
library("readxl")
library("ggalluvial")
library("foreach")
library("ggtext")
source("./scripts/common/ggpubr.A4.R")

## We focus on normal early stages only.

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
            ID, AC, AN, AF
        )]) -> .;
    } -> ..RE.matching.edits.edited.only.valid.genes.only.in.normal.early.stages.only.dt;
    ##
    ## 1.2. compute the (lower- and upper-bounds of) % of RE-matching-edits-edited transcripts per gene in per sample
    {
        copy(..RE.matching.edits.edited.only.valid.genes.only.in.normal.early.stages.only.dt) -> .;
        ## compute the lower bound of % RE-matching-edit-edited transcripts
        .[, list(
            percentage.of.RE.matching.edit.edited.transcripts.lower.bound=max(AF),
            percentage.of.RE.matching.edit.edited.transcripts.upper.bound=(1-prod(1-AF))
            ), list(SAMPLE, gse, stage, Gene_ID, Gene_Name)] -> .;
    } -> ..bounds.of.RE.matching.edit.edited.transcripts.per.gene.per.sample.valid.genes.only.in.normal.early.stages.only.dt
    ##
} -> bounds.of.RE.matching.edit.edited.transcripts.per.gene.per.sample.valid.genes.only.in.normal.early.stages.only.dt

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

## 3. Merge % edited and FPKM, and compute the spearman correlation between them
{

    ## 3.1. merge the two tables (within-sample style)
    {
        merge(x=bounds.of.RE.matching.edit.edited.transcripts.per.gene.per.sample.valid.genes.only.in.normal.early.stages.only.dt, y=combined.gexpr.FPKM.melt.normal.early.stage.samples.only.dt,
              by.x=c("SAMPLE", "Gene_ID"), by.y=c("SAMPLE", "Gene_ID"),
              all.x=TRUE, all.y=FALSE) -> .;
    } -> ..bounds.of.RE.matching.edit.edited.transcripts.and.FPKM.per.gene.per.sample.valid.genes.only.in.normal.early.stages.only.dt

    ## 3.2. compute the spearman correlation within samples
    {
    }
    
}


## [SAMPLE, stage, Gene, RE edit (ID, editing level)] 
## |
## SAMPLE, stage, Gene, % RE-edited transcripts of this gene (lower bound-upper bound)
## | <- [SAMPLE, Gene, FPKM]
## SAMPLE, stage, Gene, % RE-edited transcripts of this gene (lower bound-upper bound), FPKM
## |
## stage-before, stage-after, spearman(% RE-edited transcripts, FPKM | across all genes (could be classified into maternal clearance targets and others)) (should be overall negative if RE helps to decay)



## plot per-stage count of REs
{
    RE.dt  ->.;
    ##c select edit x stage only
    unique(.[, list(CHROM, POS, stage)]) ->.;
    ## get count of REs per stage
    .[, list(count.recurrent.edits=.N), list(stage)] ->.;
    ##s select stages to display
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula', 'blastocyst.late', 'trophoblast', 'ICM', 'CTB', 'STB', 'EVT', 'MTB', 'epiblast', 'hypoblast')] ->.;
    ## prettify stage
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, before.or.after.8.cell.ordered:=factor(before.or.after.8.cell, levels=unique(temp.stage.dt[, before.or.after.8.cell]))] ->.;
    ggplot(., aes(x=stage.description.ordered, y=count.recurrent.edits)) -> .;
    . + geom_bar(, stat="identity") -> .;
    . + labs(x="", y="# REs") ->.;
    . + scale_y_log10() ->.;
    . + facet_grid(~before.or.after.8.cell.ordered, scales="free_x", space="free_x") ->.;
    . + theme_pubr(base_size=10) ->.;
    . + theme(text=element_text(size=10), axis.text.x=element_text(angle=45, hjust=1)) ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/RE.count.per.stage.png",
        plot=.,
        width.r=0.45, height.r=0.2
    )
}


## plot exonic-intronic piecharts
plot.piechart.for.all.edit.types <- function(temp.dt){
    copy(temp.dt) ->.;
    ## collapse to [CHROM, POS, Gene_ID, Annotation.class, stage], factoring out [SAMPLE]
    ## equivalent to accounting for all edits that are observed in >=1 sample
    .[, list(CHROM, POS, Gene_ID, Annotation.class, stage)] ->.;
    unique(.) ->.;
    ## compute the count
    .[, list(count=.N), list(Annotation.class, stage)] ->.;
    ## focus on core stages only
    .[stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell", "8-cell")] ->.;
    ## prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## prettify Annotations
    .[, Annotation.class.description:=c('exonic.or.splicing.related'='exonic', 'purely.intronic'='intronic')[Annotation.class]] ->.;
    ggplot(., aes(x=1, y=count, fill=Annotation.class.description)) +
        geom_bar(stat='identity', position='fill') +
        coord_polar(theta="y") +
        facet_wrap(~stage.description.ordered, nrow=1) ->.;
    . +
        theme_pubr() +
        theme(axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
        labs(x="", y="", fill="") +
        scale_fill_manual(values=c("#F8766D", "grey50")) ->.;
    .
}

{
    observed.edit.all.edit.types.ggplot <- plot.piechart.for.all.edit.types(observed.edits.valid.genes.only.dt)
    RE.all.edit.types.ggplot <- plot.piechart.for.all.edit.types(RE.valid.genes.only.dt)
    combined.all.edit.types.ggplot <- ggarrange(
        plotlist=list(
            observed.edit.all.edit.types.ggplot + theme(legend.position="none") + ggtitle("All edits, with <= 50% <span style='color:#f8766d'>exonic</span> in most stages") + theme(plot.title = element_markdown()),
            RE.all.edit.types.ggplot + theme(legend.position="none") + ggtitle("REs, with ~ 75% <span style='color:#f8766d'>exonic</span> for all stages") + theme(plot.title=element_markdown(), plot.margin=margin(0,0,0,0.3,unit="cm")) 
        ), nrow=2, align="hv")
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/all.edit.types.piechart.png",
        plot=combined.all.edit.types.ggplot,
        width.r=0.9, height.r=0.3
    )
}


## plot 3'-UTR barplots

{
    copy(RE.valid.genes.only.dt) ->.;
    ## collapse the Annotation.corrected, grouping all non-3'-UTR exonic edits into one
    .[, Annotation.corrected.collapsed:=c("3_prime_UTR_variant"="3'-UTR", "intron_variant"="intronic")[Annotation.corrected]] ->.;
    .[is.na(Annotation.corrected.collapsed)==TRUE, Annotation.corrected.collapsed:="other"] ->.;
    ## keep exonic variant only
    .[Annotation.class=="exonic.or.splicing.related"] ->.;
    ## collapse to [CHROM, POS, Gene_ID, Annotation.corrected.collapsed, stage], factoring out [SAMPLE]
    ## equivalent to accounting for all edits that are observed in >=1 sample
    .[, list(CHROM, POS, Gene_ID, Annotation.corrected.collapsed, stage)] ->.;
    unique(.) ->.;
    ## compute the count
    .[, list(count=.N), list(Annotation.corrected.collapsed, stage)] ->.;
    ## focus on core stages only
    .[stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell", "8-cell")] ->.;
    ## prettify Annotation.corrected.collapsed
    .[, Annotation.corrected.collapsed.reordered:=factor(Annotation.corrected.collapsed, levels=c("other", "3'-UTR"))] ->.;
    ## prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;    
    ggplot(., aes(x=stage.description.ordered, y=count, fill=Annotation.corrected.collapsed.reordered)) -> .;
    . + geom_bar(stat='identity', position='fill') ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="right")  ->.;
    . + labs(x="", y="Ratio", fill="Type of\nexonic REs") -> .;
    . + scale_fill_manual(values=c("grey50", "#F8766D")) ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/3UTR.ratio.barplot.png",
        plot=., 
        width.r=0.45, height.r=0.2
    )
}


## additional: plot 3'-UTR barplots for all edits

{
    copy(observed.edits.valid.genes.only.dt) ->.;
    ## collapse the Annotation.corrected, grouping all non-3'-UTR exonic edits into one
    .[, Annotation.corrected.collapsed:=c("3_prime_UTR_variant"="3'-UTR", "intron_variant"="intronic")[Annotation.corrected]] ->.;
    .[is.na(Annotation.corrected.collapsed)==TRUE, Annotation.corrected.collapsed:="other"] ->.;
    ## keep exonic variant only
    .[Annotation.class=="exonic.or.splicing.related"] ->.;
    ## collapse to [CHROM, POS, Gene_ID, Annotation.corrected.collapsed, stage], factoring out [SAMPLE]
    ## equivalent to accounting for all edits that are observed in >=1 sample
    .[, list(CHROM, POS, Gene_ID, Annotation.corrected.collapsed, stage)] ->.;
    unique(.) ->.;
    ## compute the count
    .[, list(count=.N), list(Annotation.corrected.collapsed, stage)] ->.;
    ## focus on core stages only
    .[stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell", "8-cell")] ->.;
    ## prettify Annotation.corrected.collapsed
    .[, Annotation.corrected.collapsed.reordered:=factor(Annotation.corrected.collapsed, levels=c("other", "3'-UTR"))] ->.;
    ## prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;    
    ggplot(., aes(x=stage.description.ordered, y=count, fill=Annotation.corrected.collapsed.reordered)) -> .;
    . + geom_bar(stat='identity', position='fill') ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="right")  ->.;
    . + labs(x="", y="Ratio", fill="Type of\nexonic REs") -> .;
    . + scale_fill_manual(values=c("grey50", "#F8766D")) ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/3UTR.ratio.barplot.for.observed.edits.png",
        plot=., 
        width.r=0.45, height.r=0.2
    )
}



## plot sankey plots for pairs of early stages
## for logic consistency, here we plot REs on valid genes only


RE.CJ.on.valid.genes.only.dt <- RE.CJ.dt[CHROM.and.POS %in% RE.valid.genes.only.dt[, paste(sep="", CHROM, "_", POS)]]

temp.intermediate.dt <- {
    RE.CJ.on.valid.genes.only.dt ->.;
    ## prettify recurrence type
    .[, recurrence.type.renamed:=c("not detected"="not detected", "detected but not recurrent"="detected, not RE", "recurrent"="detected, RE")[recurrence.type]]
    .[, recurrence.type.renamed.ordered:=factor(recurrence.type.renamed, levels=c("not detected", "detected, not RE", "detected, RE"))] ->.;
    ## prettify stages
    .[, stage:=sub(pattern="@.*", replacement="", x=group)]
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## split into different stage pairs
    foreach(
        temp.stage.set=list(
            c('oocyte.GV', 'oocyte.MII'),
            c('oocyte.MII', 'zygote'),
            c('zygote', '2-cell'),
            c('2-cell', '4-cell'),
            c('4-cell', '8-cell')
        ),
        temp.part=c(
            'GV -> MII',
            'MII -> zygote',
            'zygote -> 2-cell',
            '2-cell -> 4-cell',
            '4-cell -> 8-cell'
        )) %do% {
            .[stage %in% temp.stage.set] -> temp.dt
            temp.dt[CHROM.and.POS %in% temp.dt[recurrence.type.renamed == "detected, RE", CHROM.and.POS]] -> temp.dt
            temp.dt[, part:=temp.part]
        } ->.; rbindlist(.) ->.;
    .[, part.ordered:=factor(part, levels=unique(part))]
    .
}

{
    temp.intermediate.dt -> .
    ggplot(., aes(x = stage.description.ordered, stratum = recurrence.type.renamed.ordered, alluvium = CHROM.and.POS, fill=recurrence.type.renamed.ordered)) -> .;
    . + geom_flow(stat = "flow", color = "darkgray", na.rm=FALSE) ->.;
    . + geom_stratum() ->.;
    . + facet_wrap(~ part.ordered, scales="free", nrow=1) ->.;
    . + theme_pubr(base_size=11) -> .;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) ->.;
    . + guides(fill=guide_legend(nrow=1)) ->.;
    . + scale_fill_brewer(type = "qual", palette = "Set2") ->.;
    . + labs(x="", y="# editing sites\non protein-coding genes", fill="") -> .;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/RE.sankey.png",
        plot=., 
        width.r=0.9, height.r=0.25
    )
}
