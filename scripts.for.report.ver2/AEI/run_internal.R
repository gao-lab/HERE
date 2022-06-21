library("data.table")
library("glue")
library("readxl")
library("ggtext")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/AEI"

AEI.with.ADAR.FPKM.and.sample.info.dt <- fread(glue("{output.directory}/AEI.with.ADAR.FPKM.and.sample.info.dt.gz"))

## 1. Plot AEI boxplot
{
    ##
    copy(AEI.with.ADAR.FPKM.and.sample.info.dt) -> .;
    ## plot AEI boxplot
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] -> .;
    . -> ..to.plot.dt
    ggplot(data.table(..to.plot.dt, facet="control index is not 0"), aes(x=stage.description.ordered, y=SNR)) ->.;
    . + geom_boxplot() ->.;
    . + geom_bar(data=..to.plot.dt[, list(log10.count=log10(.N), facet="control index is 0"), list(stage.description.ordered)], mapping=aes(x=stage.description.ordered, y=log10.count), stat="identity") -> .;
    . + coord_flip() -> .;
    . + facet_grid(~facet, scales="free_x") -> .;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) ->.;
    . + labs(x="", y="(boxplot) Alu editing index (AEI) = unnormalized AEI/C-to-T control index\n(barplot) log10(count of samples whose control index is 0)") ->.;
    ggsave.A4(
        filename=glue("{output.directory}/AEI.per.stage.boxplot.png"),
        plot=.,
        width.r=1, height.r=0.8
    )
    ##
}

## 2. Plot correlation with ADAR
{
    ##
    copy(AEI.with.ADAR.FPKM.and.sample.info.dt) -> .;
    melt(data=., id.vars=c("SAMPLE", "SNR", "stage", "is.normal"), measure.vars=c("ENSG00000160710.17", "ENSG00000197381.16"), variable.name="gene.id", value.name="FPKM") -> .;       
    .[is.normal==TRUE & stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell", "8-cell", "morula", "blastocyst.late", "ICM", "trophoblast", "STB", "CTB", "MTB", "EVT", "epiblast", "hypoblast", "hESC")] -> .; ## normal stages with >=10 normal samples
    .[, list(spearman.cor=cor(SNR, FPKM, method="spearman"), spearman.cor.test.pvalue=cor.test(SNR, FPKM, method="spearman", alternative="greater")$p.value), list(stage, is.normal, gene.id, gene.name=c("ENSG00000160710.17"="ADAR1", "ENSG00000197381.16"="ADAR2")[gene.id])] -> .;    
    .[, spearman.cor.test.pvalue.adjusted:=p.adjust(spearman.cor.test.pvalue, method="BH"), list(gene.name)] -> .;
    .[spearman.cor.test.pvalue.adjusted >= 0.05, spearman.cor.to.plot:=-0.6] -> .; ## just for positioning
    .[spearman.cor.test.pvalue.adjusted < 0.05, spearman.cor.to.plot:=spearman.cor] -> .; 
    ##
    ## plot AEI boxplot
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] -> .;
    . -> ..to.plot.dt
    ##
    ggplot(..to.plot.dt, aes(x=spearman.cor.to.plot, fill=paste(sep="", "<i>", gene.name, "</i>"))) ->.;
    . + geom_histogram(position="dodge", binwidth=0.3) -> .;    
    . + theme_pubr() ->.;
    . + scale_x_continuous(breaks=c(0, 0.3, 0.6, 0.9), limits=c(-0.9, 0.9)) -> .;
    . + scale_fill_manual(values=c("#303030", "#EBB130")) -> .;
    . + labs(x="<span> </span><span> </span><span> </span><span> </span>Spearman correlation<br/>between AEI and <i>ADAR</i> FPKM", y="No. of normal stages\nwith >=10 normal samples", fill="") ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), legend.text=element_markdown(), axis.title.x=element_textbox_simple(width=unit(7, "cm"), hjust=-0.03)) ->.;
    . + annotate("text", x=-0.6, y=-1, label="NS") -> .;
    . + coord_cartesian(ylim=c(0, 18), clip="off", expand=FALSE) -> .;
    ggsave.A4(
        filename=glue("{output.directory}/AEI.ADAR.correlation.histogram.png"),
        plot=.,
        width.r=0.37, height.r=0.4
    )
    ##
}
