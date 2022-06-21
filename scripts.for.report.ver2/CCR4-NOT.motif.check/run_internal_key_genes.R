library("data.table")
library("readxl")
library("magrittr")
library("glue")
library("ggpubr")
library("ggtext")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RNA.degradation.processes/"
dir.create(output.directory, recursive=TRUE)

fread("./scripts.for.report.ver2/CCR4-NOT.motif.check/key.genes.to.check.csv") -> key.genes.dt

## {
##     "./external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf" -> .;
##     fread(glue('cat {.} | grep -v "#" |grep -P "\tgene\t" | sed -E \'s/.*gene_id "([^;]+)";.* gene_name "([^;]+)".*/\\1\t\\2/\' '), header=FALSE, col.names=c("gene.id", "gene.name")) -> .;
##     .
## } -> gene.id.and.name.dt

observed.edits.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.observed.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

merge(x=observed.edits.valid.genes.only.dt, y=key.genes.dt, by.x="Gene_Name", by.y="gene.symbol", all=FALSE) -> observed.edits.key.genes.only.dt;

fwrite(observed.edits.key.genes.only.dt, glue("{output.directory}/observed.edits.key.genes.only.dt.gz"))

RE.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

merge(x=RE.valid.genes.only.dt, y=key.genes.dt, by.x="Gene_Name", by.y="gene.symbol", all=FALSE) -> RE.key.genes.only.dt;

fwrite(RE.key.genes.only.dt, glue("{output.directory}/RE.key.genes.only.dt.gz"))


MBS.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA.intersection/edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.summary.gene.and.edit.level.dt.gz")

merge(x=RE.key.genes.only.dt, y=MBS.dt, by.x=c("Gene_ID", "POS"), by.y=c("gene.id", "edit.POS"), all.x=TRUE, all.y=FALSE) -> RE.key.genes.only.with.MBS.dt

fwrite(RE.key.genes.only.with.MBS.dt, glue("{output.directory}/RE.key.genes.only.with.MBS.dt.gz"))




{
    
    RE.key.genes.only.with.MBS.dt -> .;
    ## pick edited samples
    .[is.na(AC) == FALSE] -> .;
    ## pick normal early stages
    .[ (is.normal==TRUE) & (stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell", "8-cell"))] -> .;
    ## add miRNA.type
    .[, miRNA.type:=c("TRUE;FALSE"="MBS-gaining", "FALSE;FALSE"="MBS-neutral", "FALSE;TRUE"="MBS-losing")[paste(sep="", gains.miRNA.sites, ";", loses.miRNA.sites)] ] -> .;
    ## write BED file for UCSC visualization
    fwrite(unique(.[Gene_Name=="EXOSC6", list(CHROM=CHROM, start=POS-1, end=POS-1+1, name=miRNA.type, score=1000, strand=".", thickStart=POS-1, thickEnd=POS-1+1, itemRgb=c("MBS-gaining"="237,125,49", "MBS-losing"="112,173,71", "MBS-neutral"="68,114,196")[miRNA.type])]), glue("{output.directory}/EXOSC6.RE.bed"), sep="\t", row.names=FALSE, col.names=FALSE)
    ##= prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## pick relevant columns
    .[, list(ID, AF, Gene_Name, process, stage.description.ordered, miRNA.type)] -> .;
    . -> ..to.plot.dt

    ## run ggplot
    ggplot(..to.plot.dt[Gene_Name == "EXOSC6"], aes(x=ID, y=AF, fill=miRNA.type)) -> .;
    . + geom_boxplot() ->.;
    . + coord_flip() -> .;
    . + facet_grid(~stage.description.ordered, scales="free_y", space="free_y") -> .;
    . + labs(x="", y="editing level", fill="") -> .;
    . + theme_pubr(base_size=10) ->.;
    . + scale_y_continuous(limits=c(0.1, 1), breaks=seq(0.1, 1, 0.3)) -> .;
    . + scale_fill_manual(values=c("#ED7D31", "#70AD47", "#4472C4")) -> .;
    . + ggtitle(expression(paste("REE editing profile on ", italic("EXOSC6")))) -> .;
    ggsave.A4(
        filename=glue("{output.directory}/REE.on.EXOSC6.in.normal.early.stages.AF.boxplot.png"),
        plot=.,
        width.r=0.9, height.r=0.4
    )

    ggplot(..to.plot.dt[Gene_Name == "CNOT6"], aes(x=ID, y=AF, fill=miRNA.type)) -> .;
    . + geom_boxplot() ->.;
    . + coord_flip() -> .;
    . + facet_grid(~stage.description.ordered, scales="free_y", space="free_y") -> .;
    . + labs(x="", y="editing level", fill="") -> .;
    . + theme_pubr(base_size=10) ->.;
    . + scale_y_continuous(limits=c(0.1, 1), breaks=seq(0.1, 1, 0.3)) ->.;
    . + scale_fill_manual(values=c("#4472C4")) -> .;
    . + ggtitle(expression(paste("REE editing profile on ", italic("CNOT6")))) -> .;
    ggsave.A4(
        filename=glue("{output.directory}/REE.on.CNOT6.in.normal.early.stages.AF.boxplot.png"),
        plot=.,
        width.r=0.5, height.r=0.17
    )

    
}


## Use the following setting to visualize the UCSC track
## Track file to upload: ./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RNA.degradation.processes/EXOSC6.RE.bed
## assembly: hg38
## region: chr16:70,247,600-70,250,600
## image width: 780 pixels
## text size: 12
## font: Helvetica
## Custom track: itemRgb='On', pack
## Other tracks: Base position, dense; All GENCODE (All GENCODE V32), full; RepeatMasker, full
## Track description: MBS-related REEs
