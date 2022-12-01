library("data.table")
library("readxl")
library("ggalluvial")
library("foreach")
library("ggtext")
source("./scripts/common/ggpubr.A4.R")

RE.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.dt.txt.gz")

RE.CJ.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.site.recurrence.comparison.CJ.dt.txt.gz")

observed.edits.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.observed.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

RE.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")


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
    . -> to.plot.dt;
    fwrite(to.plot.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/source.data.for.REE.count.per.stage.png.csv.gz")
    
    to.plot.dt -> .;
    ggplot(., aes(x=stage.description.ordered, y=count.recurrent.edits)) -> .;
    . + geom_bar(, stat="identity") -> .;
    . + labs(x="", y="# REEs") ->.;
    . + scale_y_log10() ->.;
    . + facet_grid(~before.or.after.8.cell.ordered, scales="free_x", space="free_x") ->.;
    . + theme_pubr(base_size=10) ->.;
    . + theme(text=element_text(size=10), axis.text.x=element_text(angle=45, hjust=1)) ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/REE.count.per.stage.png",
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
    . -> temp.to.plot.dt
    temp.to.plot.dt -> .;
    ggplot(., aes(x=1, y=count, fill=Annotation.class.description)) +
        geom_bar(stat='identity', position='fill') +
        coord_polar(theta="y") +
        facet_wrap(~stage.description.ordered, nrow=1) ->.;
    . +
        theme_pubr() +
        theme(axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
        labs(x="", y="", fill="") +
        scale_fill_manual(values=c("#F8766D", "grey50")) ->.;
    . -> temp.piechart.ggplot
    return(list(temp.to.plot.dt=temp.to.plot.dt, temp.piechart.ggplot=temp.piechart.ggplot))
}

{
    observed.edit.all.edit.types.to.plot.dt.and.ggplot.list <- plot.piechart.for.all.edit.types(observed.edits.valid.genes.only.dt)
    RE.all.edit.types.to.plot.dt.and.ggplot.list <- plot.piechart.for.all.edit.types(RE.valid.genes.only.dt)
    observed.edit.all.edit.types.to.plot.dt <- observed.edit.all.edit.types.to.plot.dt.and.ggplot.list[["temp.to.plot.dt"]]
    RE.all.edit.types.to.plot.dt <- RE.all.edit.types.to.plot.dt.and.ggplot.list[["temp.to.plot.dt"]]
    
    rbindlist(list(
        data.table(observed.edit.all.edit.types.to.plot.dt, subpanel="observed"),
        data.table(RE.all.edit.types.to.plot.dt, subpanel="REE")
    ), use.names=TRUE) -> temp.combined.dt
    fwrite(temp.combined.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/source.data.for.REE.all.edit.types.piechart.png.csv.gz")
    
    ##
    observed.edit.all.edit.types.ggplot <- observed.edit.all.edit.types.to.plot.dt.and.ggplot.list[["temp.piechart.ggplot"]]
    RE.all.edit.types.ggplot <- RE.all.edit.types.to.plot.dt.and.ggplot.list[["temp.piechart.ggplot"]]
    combined.all.edit.types.ggplot <- ggarrange(
        plotlist=list(
            observed.edit.all.edit.types.ggplot + theme(legend.position="none") + ggtitle("All edits, with <= 50% <span style='color:#f8766d'>exonic</span> in most stages") + theme(plot.title = element_markdown()),
            RE.all.edit.types.ggplot + theme(legend.position="none") + ggtitle("REEs, with ~ 75% <span style='color:#f8766d'>exonic</span> for all stages") + theme(plot.title=element_markdown(), plot.margin=margin(0,0,0,0.3,unit="cm")) 
        ), nrow=2, align="hv")
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/REE.all.edit.types.piechart.png",
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
    . -> to.plot.dt
    fwrite(to.plot.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/source.data.for.REE.3UTR.ratio.barplot.png.csv.gz")

    to.plot.dt -> .;
    ggplot(., aes(x=stage.description.ordered, y=count, fill=Annotation.corrected.collapsed.reordered)) -> .;
    . + geom_bar(stat='identity', position='fill') ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="right")  ->.;
    . + labs(x="", y="Ratio", fill="Type of\nexonic REEs") -> .;
    . + scale_fill_manual(values=c("grey50", "#F8766D")) ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/REE.3UTR.ratio.barplot.png",
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
    . + labs(x="", y="Ratio", fill="Type of\nexonic edits") -> .;
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
    .[, recurrence.type.renamed:=c("not detected"="not detected", "detected but not recurrent"="detected, not REE", "recurrent"="detected, REE")[recurrence.type]]
    .[, recurrence.type.renamed.ordered:=factor(recurrence.type.renamed, levels=c("not detected", "detected, not REE", "detected, REE"))] ->.;
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
            temp.dt[CHROM.and.POS %in% temp.dt[recurrence.type.renamed == "detected, REE", CHROM.and.POS]] -> temp.dt
            temp.dt[, part:=temp.part]
        } ->.; rbindlist(.) ->.;
    .[, part.ordered:=factor(part, levels=unique(part))]
    .
}

fwrite(temp.intermediate.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/source.data.for.REE.sankey.png.csv.gz")

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
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE/REE.sankey.png",
        plot=., 
        width.r=0.9, height.r=0.25
    )
}
