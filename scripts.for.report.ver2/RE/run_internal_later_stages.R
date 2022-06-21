library("data.table")
library("readxl")
library("ggalluvial")
library("foreach")
library("ggtext")
library("glue")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.later.stages/"
dir.create(output.directory, recursive=TRUE)

RE.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.dt.txt.gz")

RE.CJ.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.site.recurrence.comparison.CJ.dt.txt.gz")

observed.edits.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.observed.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

RE.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")



## plot exonic-intronic piecharts (for later stages)
plot.piechart.for.all.edit.types <- function(temp.dt){
    copy(temp.dt) ->.;
    ## collapse to [CHROM, POS, Gene_ID, Annotation.class, stage], factoring out [SAMPLE]
    ## equivalent to accounting for all edits that are observed in >=1 sample
    .[, list(CHROM, POS, Gene_ID, Annotation.class, stage)] ->.;
    unique(.) ->.;
    ## compute the count
    .[, list(count=.N), list(Annotation.class, stage)] ->.;
    .[, total.count:=sum(count), stage] -> .;
    ## focus on core stages only
    .[stage %in% c("morula", "blastocyst.late", "ICM", "CTB", "STB", "EVT", "MTB", "epiblast", "hypoblast")] ->.;
    ## prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## prettify Annotations
    .[, Annotation.class.description:=c('exonic.or.splicing.related'='exonic', 'purely.intronic'='intronic')[Annotation.class]] ->.;
    ## prettify count label
    .[, `:=`(count.label=paste(sep="", {.SD[Annotation.class.description=="exonic", count] -> temp.exonic.count; if (length(temp.exonic.count)==0){temp.exonic.count <- 0} ; temp.exonic.count}, "/", .SD[1, total.count])), stage] -> .;
    ##
    . -> ..to.plot.dt  
    ggplot(..to.plot.dt, aes(x=1, y=count, fill=Annotation.class.description)) +
        geom_bar(stat='identity', position='fill') +
        coord_polar(theta="y") +
        facet_wrap(~stage.description.ordered + count.label, nrow=1) ->.;
    . +
        theme_pubr(base_size=10) +
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
            observed.edit.all.edit.types.ggplot + theme(legend.position="none") + ggtitle("All edits (with <span style='color:#f8766d'>exonic</span> shown in <span style='color:#f8766d'>red</span>)") + theme(plot.title = element_markdown()),
            RE.all.edit.types.ggplot + theme(legend.position="none") + ggtitle("REEs (with <span style='color:#f8766d'>exonic</span> shown in <span style='color:#f8766d'>red</span>)") + theme(plot.title=element_markdown(), plot.margin=margin(0,0,0,0.3,unit="cm")) 
        ), nrow=2, align="hv")
    ggsave.A4(
        filename=glue("{output.directory}/REE.all.edit.types.piechart.later.stages.png"),
        plot=combined.all.edit.types.ggplot,
        width.r=0.97, height.r=0.3
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
    .[stage %in% c("morula", "blastocyst.late", "ICM", "CTB", "STB", "EVT", "MTB", "epiblast", "hypoblast")] -> .;
    ## prettify Annotation.corrected.collapsed
    .[, Annotation.corrected.collapsed.reordered:=factor(Annotation.corrected.collapsed, levels=c("other", "3'-UTR"))] ->.;
    ## prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;    
    ggplot(., aes(x=stage.description.ordered, y=count, fill=Annotation.corrected.collapsed.reordered)) -> .;
    . + geom_bar(stat='identity', position='fill') ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="right")  ->.;
    . + labs(x="", y="Ratio", fill="Type of\nexonic REEs") -> .;
    . + scale_fill_manual(values=c("grey50", "#F8766D")) ->.;
    ggsave.A4(
        filename=glue("{output.directory}/REE.3UTR.ratio.barplot.later.stages.png"),
        plot=., 
        width.r=0.45, height.r=0.2
    )
}
