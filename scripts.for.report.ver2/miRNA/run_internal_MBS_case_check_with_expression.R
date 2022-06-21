library("data.table")
library("glue")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA.case.study/"
dir.create(output.directory, recursive=TRUE)

## 1. Get editing level of REs of SUV39H per sample in early stage
{
    
    ## 1.1. Get REs on SUV39H2 with MBS annotation
    {
        fread(
            "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/case.study.bed", header=FALSE,
            col.names=c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb")
        ) -> .;
        .[, list(CHROM=chrom, POS=chromStart + 1, MBS.type=name)] -> .;
        .
    } -> ..temp.SUV39H2.RE.MBS.info.dt
    
    ##
    ## 1.2. Get editing level for each RE edited events on valid genes in normal early stages
    {
        fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz") -> .;
        ## keep edited events only
        .[is.na(AC) == FALSE] -> .;
        ## keep normal early stages only
        .[is.normal == TRUE & stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell')] -> .;
        ## merge with SVU39H2 RE MBS info
        merge(x=., y=..temp.SUV39H2.RE.MBS.info.dt,
              by.x=c("CHROM", "POS"), by.y=c("CHROM", "POS"),
              all.x=FALSE, all.y=FALSE) -> .;
        .[, list(ID, CHROM, POS, MBS.type, SAMPLE, gse, stage, AF)] -> .;
    } -> ..temp.SUV39H2.RE.MBS.and.AF.per.sample.info.dt
    fwrite(..temp.SUV39H2.RE.MBS.and.AF.per.sample.info.dt, glue("{output.directory}/temp.SUV39H2.RE.MBS.and.AF.per.sample.info.dt.gz"))
    ..temp.SUV39H2.RE.MBS.and.AF.per.sample.info.dt
} -> temp.SUV39H2.RE.MBS.and.AF.per.sample.info.dt

fwrite(unique(temp.SUV39H2.RE.MBS.and.AF.per.sample.info.dt[, list(CHROM=CHROM, start=POS-1, end=POS-1+1, name=MBS.type, score=1000, strand=".", thickStart=POS-1, thickEnd=POS-1+1, itemRgb=c("MBS-gaining"="237,125,49", "MBS-losing"="112,173,71", "MBS-neutral"="68,114,196")[MBS.type])]), glue("{output.directory}/case.study.key.stage.only.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Visualization by UCSC:
## region: chr10:14,903,970-14,904,120
## tracks displayed: the "My track" (full), Base position (dense), and ALL GENCODE v32 (full)
## image width: 720 pixels; label area width: 20 characters; text size: 14; font: Helvetica; style: Normal

## 2. Get FPKM of SUV39H2 per sample in early stage 
{
    ##
    ## 2.1. get sample info of all normal oocyte GV stage samples (note that the REs on SUV39H2 in normal samples are all in oocyte GV)
    {
        fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt") -> .;
        ## keep normal early stage samples
        .[is.normal==TRUE & stage == 'oocyte.GV'] -> .;
    } -> ..normal.early.stage.samples.info.dt
    ##
    ## 2.2. get SUV39H2 FPKM for all normal early stage samples (only samples in the 2,071 set are considered; the rest are discarded due to insufficient read length (e.g, GSM1160130, is 8-cell, but has read length of 49*2 < 75*2)
    {
        read.table("result/BS06_1__get_expression_level/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/combined.gexpr.FPKM.matrix.txt") -> .;
        ## keep matched samples only
        .[, intersect(colnames(.), ..normal.early.stage.samples.info.dt[, gsm])] -> .;
        ## transform to 3-way table
        data.table(Gene_ID=rownames(.), .) -> .;
        melt(., id.vars="Gene_ID", variable.name="SAMPLE", value.name="FPKM") -> .;
        ## pick SUV39H2 expression only
        .[Gene_ID=="ENSG00000152455.16"] -> .;
    } -> ..SUV39H2.FPKM.melt.normal.early.stage.samples.only.dt;
    fwrite(..SUV39H2.FPKM.melt.normal.early.stage.samples.only.dt, glue("{output.directory}/SUV39H2.FPKM.melt.normal.early.stage.samples.only.dt.gz"))
    ##
    ..SUV39H2.FPKM.melt.normal.early.stage.samples.only.dt
} -> SUV39H2.FPKM.melt.normal.early.stage.samples.only.dt


## 3. merge the editing level table with the FPKM table
{
    merge(
        x=temp.SUV39H2.RE.MBS.and.AF.per.sample.info.dt,
        y=SUV39H2.FPKM.melt.normal.early.stage.samples.only.dt,
        by.x="SAMPLE", by.y="SAMPLE",
        all.x=TRUE, all.y=FALSE ## discard samples with no such RE detected
    ) -> ..SUV39H2.RE.MBS.and.AF.and.FPKM.per.sample.normal.early.stage.samples.only.dt;
} -> SUV39H2.RE.MBS.and.AF.and.FPKM.per.sample.normal.early.stage.samples.only.dt

## 4. plot the trend
{
    
    copy(SUV39H2.RE.MBS.and.AF.and.FPKM.per.sample.normal.early.stage.samples.only.dt) -> .;
    .[, spearman.cor:=cor(AF, FPKM, method="spearman"), list(ID)] -> .;
    . -> ..to.plot.dt;
    ggplot(..to.plot.dt, aes(x=AF, y=FPKM, color=MBS.type)) -> .;
    . + geom_point(size = 3) -> .;
    . + geom_smooth(method=lm, color="black") -> .;
    . + geom_text(data=unique(..to.plot.dt[, list(ID, MBS.type, spearman.cor, AF=0.18, FPKM=12)])[ID=='chr10_14904041_A_G', `:=`(AF=0.23, FPKM=23)][ID=="chr10_14904104_A_G", `:=`(AF=0.16, FPKM=21)], mapping=aes(label=glue("Spearman's rho:\n{round(spearman.cor, 4)}")), size=3, color="black") -> .;
    ##. + scale_x_continuous(limits=c(0.1, 1), breaks=c(0.1, 0.4, 0.7, 1)) -> .;
    ##. + scale_y_log10(limits=c(1e-3, 1e4), breaks=c(1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4), labels=function(x) {sub(pattern="(\\.0+|0+)$", replacement="", x=format(x, scientific=FALSE))}) -> .;
    . + theme_pubr(base_size=12) -> .;
    ##. + theme(axis.text.x=element_text(angle=45, hjust=1)) -> .;
    . + facet_wrap(~ID, scales="free", nrow=1) -> .;
    . + labs(x="editing level", y="FPKM", color="") -> .;
    ggsave.A4(filename=glue("{output.directory}/SUV39H2.AF.vs.FPKM.scatterplot.oocyte.GV.png"), plot=., width.r=0.9, height.r=0.3)

}
