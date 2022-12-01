library("data.table")
library("readxl")
library("ggalluvial")
library("foreach")
library("magrittr")
source("./scripts/common/ggpubr.A4.R")

RE.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

total.sample.count.for.normal.stages.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.sample.count.for.normal.stages.dt.csv")
total.sample.count.for.normal.stages.mapping.vector <- {
    total.sample.count.for.normal.stages.dt -> .;
    temp.vector <- .[, total.sample.count];
    names(temp.vector) <- .[, stage]
    temp.vector
}

tables.to.save.list <- list()

## generate table of RE-targeted genes and plot counts per stage
{

    RE.valid.genes.only.dt -> .;
    ##s pick edited records only
    .[is.na(AC)==FALSE] ->.;
    ##c collapse Annotation.class 
    .[, list(Annotation.class.pasted=paste(collapse=";", Annotation.class %>% sort %>% unique)), list(Gene_ID, Gene_Name, SAMPLE, gse, stage)] ->.;
    .[, Annotation.class.pasted.description:=c('exonic.or.splicing.related'='exonic', 'purely.intronic'='intronic', 'exonic.or.splicing.related;purely.intronic'='mixed')[Annotation.class.pasted]] ->.;
    ##= compute % occurrence of each gene x Annotation.class.pasted.
    .[, total.sample.count:=total.sample.count.for.normal.stages.mapping.vector[stage]]
    .[, gene.occurrence:=.N, list(Gene_ID, Gene_Name, stage)]
    .[, gene.occurrence.pct:=gene.occurrence/total.sample.count]
    ##c pick those with gene occurrence pct >= 80% samples. Note that the resulting `dt` conforms to alluvial form (i.e., <=1 record for each gene x stage)
    .[gene.occurrence.pct>=0.8] ->.;
    ##c collapse to each stage
    .[, list(Annotation.class.pasted.description.mixed=paste(collapse=";", Annotation.class.pasted.description %>% sort %>% unique)), list(stage, Gene_ID, Gene_Name)] ->.;
    ##= reduce Annotation to exon/intron/mixed
    .[, Annotation.class.pasted.description.mixed.ordered:=factor(c("exonic"="always exonic\n(mostly 3'-UTR)", "intronic"="always intronic", "mixed"="mixed", "exonic;intronic"="mixed", "exonic;mixed"="mixed", "intronic;mixed"="mixed", "exonic;intronic;mixed"="mixed")[Annotation.class.pasted.description.mixed], levels=c("always exonic\n(mostly 3'-UTR)", "always intronic", "mixed"))] 
    ##s only listing those 'normal' stages with sample size >= 10
    ## .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula', 'blastocyst.late', 'trophoblast', 'ICM', 'CTB', 'STB', 'MTB', 'EVT', 'epiblast', 'hypoblast')] ->.;
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula')] ->.;
    tables.to.save.list[["RE.targeted.genes.dt"]] <- copy(.)
    ## prettify stage
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, before.or.after.8.cell.ordered:=factor(before.or.after.8.cell, levels=unique(temp.stage.dt[, before.or.after.8.cell]))] ->.;
    . -> to.plot.dt
    
    fwrite(to.plot.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.gene/source.data.for.REE.gene.count.per.stage.png.csv.gz")

    ## start plotting
    to.plot.dt -> .;    
    ggplot(., aes(x = stage.description.ordered, fill=Annotation.class.pasted.description.mixed.ordered)) -> .;
    . + geom_bar(position="stack") ->.;
    . + scale_fill_manual(values=c("#FFBDC8", "#A0C5E3", "#C99CCA", "gray50")) ->.;
    ## . + facet_wrap(~before.or.after.8.cell.ordered, scales="free") ->.;
    . + labs(x="", y="# REE-targeted genes", fill="within-stage\nediting style") ->.;
    . + theme_pubr(base_size=10) ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="right") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.gene/REE.gene.count.per.stage.png",
        plot=.,
        width.r=0.45, height.r=0.2
    )    
}
fwrite(tables.to.save.list[["RE.targeted.genes.dt"]], "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.gene/RE.targeted.genes.dt.csv.gz")

## generate sankey plot
{
    
    copy(tables.to.save.list[["RE.targeted.genes.dt"]]) ->.;
    ##s pick stages of interest
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell')] ->.;
    ##a add missing gene-stage combination by CJ
    setkey(., stage, Gene_ID)
    .[CJ(stage, Gene_ID, unique=TRUE)] ->.;
    ## split into different stage pairs
    .[is.na(Annotation.class.pasted.description.mixed) == TRUE, gene.recurrently.edited:=FALSE]
    .[is.na(Annotation.class.pasted.description.mixed) == FALSE, gene.recurrently.edited:=TRUE]
    foreach(
        temp.stage.set=list(
            c('oocyte.GV', 'oocyte.MII'),
            c('oocyte.MII', 'zygote'),
            c('zygote', '2-cell'),
            c('2-cell', '4-cell')
        ),
        temp.part=c(
            'GV -> MII',
            'MII -> zygote',
            'zygote -> 2-cell',
            '2-cell -> 4-cell'
        )) %do% {
            .[stage %in% temp.stage.set] -> temp.dt
            temp.dt[Gene_ID %in% temp.dt[gene.recurrently.edited == TRUE, Gene_ID]] -> temp.dt
            temp.dt[is.na(Annotation.class.pasted.description.mixed) == TRUE, Annotation.class.pasted.description.mixed:="not REE-targeted"]
            temp.dt[, part:=temp.part]
        } ->.; rbindlist(.) ->.;
    .[, part.ordered:=factor(part, levels=unique(part))]
    ## re-prettify Annotation
    .[, Annotation.class.pasted.description.mixed.ordered:=factor(c("exonic"="always exonic\n(mostly 3'-UTR)", "intronic"="always intronic", "mixed"="mixed", "exonic;intronic"="mixed", "exonic;mixed"="mixed", "intronic;mixed"="mixed", "exonic;intronic;mixed"="mixed", "not REE-targeted"="not REE-targeted")[Annotation.class.pasted.description.mixed], levels=c("always exonic\n(mostly 3'-UTR)", "always intronic", "mixed", "not REE-targeted"))] 
    ## prettify stage
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, before.or.after.8.cell.ordered:=factor(before.or.after.8.cell, levels=unique(temp.stage.dt[, before.or.after.8.cell]))] ->.;
    . -> to.plot.dt

    fwrite(to.plot.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.gene/source.data.for.REE.gene.sankey.png.csv.gz")

    
    ## make the plot
    to.plot.dt -> .;
    ggplot(., aes(x = stage.description.ordered, stratum = Annotation.class.pasted.description.mixed.ordered, alluvium = Gene_ID, fill=Annotation.class.pasted.description.mixed.ordered)) ->.;
    . + geom_flow(stat = "flow", color = "darkgray", na.rm=FALSE) + geom_stratum() ->.;
    . + facet_wrap(~ part.ordered, scales="free", nrow=3) ->.;
    . + labs(x="", y="# genes", fill="within-stage\nediting style") ->.;
    . + scale_fill_manual(values=c("lightgreen", "#FFBDC8", "#A0C5E3", "#C99CCA", "gray50")) ->.;
    . + theme_pubr(base_size=10) ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) ->.;
    . + guides(fill=guide_legend(nrow=2)) ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.gene/REE.gene.sankey.png",
        plot=.,
        width.r=0.45, height.r=0.4
    )    
}
