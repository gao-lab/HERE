library("data.table")
library("glue")
library("readxl")
library("ggtext")
library(statpsych)
source("./scripts/common/ggpubr.A4.R")

output.directory <- "report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/AEI"

## 1. get AEI

{
    "./result/BS80_1__compute_AEI/210215-sixth-dataset/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/combined.Editing.file.dt.gz" -> .;
    fread(.) -> .;
    . -> AEI.dt
}

## 2. get ADAR expression
{
    ##
    "result/BS06_1__get_expression_level/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/combined.gexpr.FPKM.matrix.txt" -> ..TEMP.FILENAME;
    fread(..TEMP.FILENAME, header=FALSE, skip=1) -> .
    setnames(., c("gene.id", scan(..TEMP.FILENAME, character(), nlines=1))) -> .;
    ## ADARB1/ADAR2: ENSG00000197381.16
    ## ADAR1: ENSG00000160710.17
    .[gene.id %in% c("ENSG00000197381.16", "ENSG00000160710.17")] -> .;
    ## melt and dcast to get sample x gene
    melt(data=., id.vars="gene.id", variable.name="SAMPLE", value.name="FPKM") -> .;
    dcast(., SAMPLE ~ gene.id, value.var="FPKM") -> .;
    ##
    .        
} -> ADAR.FPKM.dt
fwrite(ADAR.FPKM.dt, glue("{output.directory}/ADAR.FPKM.dt.gz"))

## 3. merge the two

{
    merge(
        x=AEI.dt, y=ADAR.FPKM.dt,
        by="SAMPLE",
        all.x=TRUE, all.y=FALSE
    ) -> .;
    ##
    .
} -> AEI.with.ADAR.FPKM.dt
fwrite(AEI.with.ADAR.FPKM.dt, glue("{output.directory}/AEI.with.ADAR.FPKM.dt.gz"))

## 4. add sample info

{
    fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt") -> .;
    merge(
        x=AEI.with.ADAR.FPKM.dt, y=.,
        by.x="SAMPLE", by.y="gsm",
        all.x=TRUE, all.y=FALSE
    ) -> .;
    ##
    .
} -> AEI.with.ADAR.FPKM.and.sample.info.dt
fwrite(AEI.with.ADAR.FPKM.and.sample.info.dt, glue("{output.directory}/AEI.with.ADAR.FPKM.and.sample.info.dt.gz"))



## 1. Plot AEI boxplot
{
    ##
    copy(AEI.with.ADAR.FPKM.and.sample.info.dt) -> .;
    ## plot AEI boxplot
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] -> .;
    . -> ..to.plot.dt
    ggplot(..to.plot.dt, aes(x=stage.description.ordered, y=A2GEditingIndex)) ->.;
    . + geom_boxplot() ->.;
    . + coord_flip() -> .;
    . + theme_pubr() ->.;
    . + labs(x="", y="Alu editing index (AEI)") ->.;
    ggsave.A4(
        filename=glue("{output.directory}/official.AEI.per.stage.boxplot.png"),
        plot=.,
        width.r=1, height.r=0.8
    )
    ##
}

## 2. Plot correlation with ADAR
{
    
    copy(AEI.with.ADAR.FPKM.and.sample.info.dt) -> .;
    melt(data=., id.vars=c("SAMPLE", "A2GEditingIndex", "stage", "is.normal"), measure.vars=c("ENSG00000160710.17", "ENSG00000197381.16"), variable.name="gene.id", value.name="FPKM") -> .;       
    .[is.normal==TRUE & stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell", "8-cell", "morula", "blastocyst.late", "ICM", "trophoblast", "STB", "CTB", "MTB", "EVT", "epiblast", "hypoblast", "hESC")] -> .; ## normal stages with >=10 normal samples
    .[, list(
        spearman.n=.N,
        spearman.cor=cor(A2GEditingIndex, FPKM, method="spearman"),
        spearman.cor.test.pvalue=cor.test(A2GEditingIndex, FPKM, method="spearman", alternative="two.sided")$p.value,
        spearman.cor.test.95CI.LL=ci.spear(alpha=0.5, y=A2GEditingIndex, x=FPKM)[, "LL"],
        spearman.cor.test.95CI.UL=ci.spear(alpha=0.5, y=A2GEditingIndex, x=FPKM)[, "UL"]),
      list(stage, is.normal, gene.id, gene.name=c("ENSG00000160710.17"="ADAR1", "ENSG00000197381.16"="ADAR2")[gene.id])] -> .;    
    .[, spearman.cor.test.pvalue.adjusted:=p.adjust(spearman.cor.test.pvalue, method="BH"), list(gene.name)] -> .;
    .[spearman.cor.test.pvalue.adjusted >= 0.05, spearman.cor.to.plot:=-0.6] -> .; ## just for positioning
    .[spearman.cor.test.pvalue.adjusted < 0.05, spearman.cor.to.plot:=spearman.cor] -> .; 
    ##
    ## plot AEI boxplot
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] -> .;
    . -> ..to.plot.dt

    fwrite(..to.plot.dt, glue("{output.directory}/stat.info.for.official.AEI.ADAR.correlation.histogram.png.csv.gz"))
    
    ggplot(..to.plot.dt, aes(x=spearman.cor.to.plot, fill=paste(sep="", "<i>", gene.name, "</i>"))) ->.;
    . + geom_histogram(position="dodge", binwidth=0.3) -> .;    
    . + theme_pubr() ->.;
    ##. + scale_x_continuous(breaks=c(0, 0.3, 0.6, 0.9), limits=c(-0.9, 0.9)) -> .;
    . + scale_fill_manual(values=c("#303030", "#EBB130")) -> .;
    . + labs(x="<span> </span><span> </span><span> </span><span> </span>Spearman correlation<br/>between AEI and <i>ADAR</i> FPKM", y="No. of normal stages\nwith >=10 normal samples", fill="") ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), legend.text=element_markdown(), axis.title.x=element_textbox_simple(width=unit(7, "cm"), hjust=-0.03)) ->.;
    . + annotate("text", x=-0.6, y=-1, label="NS") -> .;
    . + coord_cartesian(ylim=c(0, 18), clip="off", expand=FALSE) -> .;
    ggsave.A4(
        filename=glue("{output.directory}/official.AEI.ADAR.correlation.histogram.png"),
        plot=.,
        width.r=0.37, height.r=0.4
    )
    ##
}
