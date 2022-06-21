library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("scales")
library("ggdendro")
library("GEOquery")
library("glue")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.RE.matching.edits.and.expression/"
dir.create(output.directory, recursive=TRUE)

## Step 1. get SAMPLE patient info
{
    ## READ GEO SOFT
    getGEO(filename="external/NCBI.GEO.SOFT/GSE95477_family.soft.gz") -> .;
    ## GET GSMList
    GSMList(.) -> .;
    ## GET GSM meta for each GSM
    rbindlist(lapply(., FUN=function(temp.GSM){
        Meta(temp.GSM) -> .a
        data.table(title=.a[["title"]], SAMPLE=.a[["geo_accession"]], developmental.stage=grep("developmental stage:", .a[["characteristics_ch1"]], value=TRUE))
    }), use.names=TRUE) -> .;
    ## ADD patient ID
    .[, patient.ID:=sub(pattern="(YNG|AMA).* - ", replacement="Patient_\\1_", x=title)] -> .;
    ## ADD stage
    .[, stage:=c("developmental stage: germinal vesicle"="oocyte.GV", "developmental stage: metaphase II"="oocyte.MII")[developmental.stage]]
    ## KEEP COLUMNS SAMPLE, patient.ID
    .[, list(SAMPLE, patient.ID, stage)] -> .;
    ## DCAST
    dcast(., patient.ID ~ stage, value.var="SAMPLE") -> .;
    setnames(., c("oocyte.GV", "oocyte.MII"), c("SAMPLE.oocyte.GV", "SAMPLE.oocyte.MII")) -> .;
    ## RETURN
    .
} -> GSE95477.patient.info.dt

## Step 2. Get expression tables for GSE95477
## START
## |<- 2.1. Determine the sample accessions to look up FPKM [..GSE95477.sample.accessions.vector]
## 2.2 Read for these samples expression table (3-way format) [..combined.gexpr.FPKM.melt.GSE95477.samples.only.dt]
## 2.2.BACKUP
## RETURN
## NOTE: do not use the edit x FPKM table, because it does not contain genes with no edits per sample
{
    ## 2.1. Determine the sample accessions to look up FPKM
    {
        c(GSE95477.patient.info.dt[, SAMPLE.oocyte.GV], GSE95477.patient.info.dt[, SAMPLE.oocyte.MII]) -> .;
        .
    } -> ..GSE95477.sample.accessions.vector
    ## 
    ## 2.2. Read for these samples expression table (3-way format)
    {
        read.table("result/BS06_1__get_expression_level/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/combined.gexpr.FPKM.matrix.txt")[, ..GSE95477.sample.accessions.vector] -> .;
        data.table(Gene_ID=rownames(.), .) -> .;
        melt(., id.vars="Gene_ID", variable.name="SAMPLE", value.name="FPKM") -> .;
    } -> ..combined.gexpr.FPKM.melt.GSE95477.samples.only.dt;
    ##
    ## 2.2.BACKUP. write current result to file
    {
        fwrite(..combined.gexpr.FPKM.melt.GSE95477.samples.only.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/combined.gexpr.FPKM.melt.GSE95477.samples.only.dt.gz")
    }
    ##
    ## RETURN
    ..combined.gexpr.FPKM.melt.GSE95477.samples.only.dt
} -> combined.gexpr.FPKM.melt.GSE95477.samples.only.dt



## Step 3. Get RE-matched edit count in GSE95477 non-outlier samples
## START
## |<- 3.1. Get RE and relevant info present in GSE95477 stages [..GSE95477.RE.and.gene.per.stage.dt]
## 3.2. Filter all GSE95477 edits for those matching the REs above [..GSE95477.RE.matched.edits.dt]
## 3.3. Compute RE-matched edit count per (Gene, SAMPLE) [..GSE95477.RE.matched.edits.count.dt]
## |<- 3.4. Get maternal gene info [..maternal.gene.info.dt]
## 3.5. Keep maternal genes only [..GSE95477.RE.matched.edits.on.maternal.gene.only.dt]
## 3.5.BACKUP.
{
    ##
    ## 3.1. Get RE and relevant info (i.e., position and their genes, and stage) present in GSE95477 stages (i.e.,oocytes GV and oocytes MII)
    ## must take unique combinations
    {
        fread(
            "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz",
            select=c("CHROM", "POS", "Gene_ID", "Gene_Name", "stage")
        )[stage %in% c('oocyte.GV', 'oocyte.MII')] -> .;
        unique(.) -> .;
        .
    } -> ..GSE95477.RE.and.gene.per.stage.dt;
    ##    
    ## 3.2. Filter all GSE95477 edits for those matching the REs above
    {
        ## read edits in all samples
        fread(
            "result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz"
        )[gse == "GSE95477"] -> .;
        ## filter for GSE95477 RE-matched edits
        merge(x=., y=..GSE95477.RE.and.gene.per.stage.dt,
              by.x=c("CHROM", "POS", "stage"), by.y=c("CHROM", "POS", "stage"),
              all.x=FALSE, all.y=FALSE, allow.cartesian=TRUE) ->.;
    } -> ..GSE95477.RE.matched.edits.dt;
    ##
    ## 3.3. Compute RE-matched edit count per (Gene, SAMPLE)
    {
        ..GSE95477.RE.matched.edits.dt[, list(RE.matched.edits.count=.N), list(Gene_ID, Gene_Name, gse, SAMPLE, stage, maternal.age, maternal.bin.ordered=factor(c("young", "old")[(as.integer(maternal.age) > 35)+1], levels=c("old", "young")))] -> .;
    } -> ..GSE95477.RE.matched.edits.count.dt;
    ##
    ## 3.4. Get maternal gene info
    {
        fread(
            "result/S42_1__annotate_embryonic_genes/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/201221-fifth-phenotype-collection/combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt.txt.gz"
        ) -> .;
        .[
            cluster %in% c("maternal.decay", "maternal.others"),
            list(Gene_ID=gene.id, Gene_Name=gene.name, cluster)
        ] -> .;
    } -> ..maternal.gene.info.dt;
    ##
    ## 3.5. Keep maternal genes only
    {
        merge(x=..GSE95477.RE.matched.edits.count.dt, y=..maternal.gene.info.dt,
          by=c("Gene_ID", "Gene_Name"),
          all.x=FALSE, all.y=FALSE) ->.;
    } -> ..GSE95477.RE.matched.edits.count.on.maternal.gene.only.dt;
    ##
    ## 3.5.BACKUP. write current result to file
    {
        fwrite(..GSE95477.RE.matched.edits.count.on.maternal.gene.only.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.RE.matched.edits.count.on.maternal.gene.only.dt.gz")
    }
    ##
    ..GSE95477.RE.matched.edits.count.on.maternal.gene.only.dt
} -> GSE95477.RE.matched.edits.count.on.maternal.gene.only.dt; 

## Step 4. Examine for each gene RE-matched edit count (oocyte GV) ~ expression (oocyte MII)
## NOTE: for each oocyte GV samples, only those genes with >=1 RE-matched edit were considered
## START
## 4.1. remove outliers in phenotype table [..GSE95477.patient.info.no.outliers.dt]
## 4.2. add RE-matched edit counts per gene to oocyte GV samples [..GSE95477.patient.info.with.oocyte.GV.RE.matched.edit.counts.dt]
## 4.3. add expression per gene to oocyte MII samples [..GSE95477.patient.info.with.oocyte.GV.RE.matched.edit.counts.and.with.oocyte.MII.FPKM.dt]
## 4.3.BACKUP.
## RETURN
{
    ##
    ## 4.1. remove outliers in phenotype table [..GSE95477.patient.info.no.outliers.dt]
    {
        GSE95477.patient.info.dt -> .;
        .[(SAMPLE.oocyte.GV %in% c("GSM2514781","GSM2514773")) == FALSE] -> .;
        .[(SAMPLE.oocyte.MII %in% c("GSM2514781","GSM2514773")) == FALSE] -> .;
        .
    } -> ..GSE95477.patient.info.no.outliers.dt
    ##
    ## 4.2. add RE-matched edit counts per gene to oocyte GV samples
    {
        merge(
            x=..GSE95477.patient.info.no.outliers.dt,
            y=GSE95477.RE.matched.edits.count.on.maternal.gene.only.dt[, list(SAMPLE.oocyte.GV=SAMPLE, maternal.bin.ordered, Gene_ID, Gene_Name, cluster, RE.matched.edits.count.oocyte.GV=RE.matched.edits.count)],
            by="SAMPLE.oocyte.GV",
            all.x=TRUE, all.y=FALSE) -> .;
    } -> ..GSE95477.patient.info.no.outliers.with.oocyte.GV.RE.matched.edit.counts.dt;
    ##
    ## 4.3. add expression per gene to oocyte MII samples
    {
        merge(
            x=..GSE95477.patient.info.no.outliers.with.oocyte.GV.RE.matched.edit.counts.dt,
            y=combined.gexpr.FPKM.melt.GSE95477.samples.only.dt[, list(SAMPLE.oocyte.MII=SAMPLE, Gene_ID, FPKM.oocyte.MII=FPKM)],
          by=c("SAMPLE.oocyte.MII", "Gene_ID"),
          all.x=TRUE, all.y=FALSE) -> .;
    } -> ..GSE95477.patient.info.no.outliers.with.oocyte.GV.RE.matched.edit.counts.and.with.oocyte.MII.FPKM.dt
    ##
    ## 4.3.BACKUP. write current result to file
    {
        fwrite(..GSE95477.patient.info.no.outliers.with.oocyte.GV.RE.matched.edit.counts.and.with.oocyte.MII.FPKM.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.patient.info.no.outliers.with.oocyte.GV.RE.matched.edit.counts.and.with.oocyte.MII.FPKM.dt.gz")
    }
    ..GSE95477.patient.info.no.outliers.with.oocyte.GV.RE.matched.edit.counts.and.with.oocyte.MII.FPKM.dt
} -> GSE95477.patient.info.no.outliers.with.oocyte.GV.RE.matched.edit.counts.and.with.oocyte.MII.FPKM.dt

## Step 5. analyze the patterns 
## GSE95477.patient.info.no.outliers.with.oocyte.GV.RE.matched.edit.counts.and.with.oocyte.MII.FPKM.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.patient.info.no.outliers.with.oocyte.GV.RE.matched.edit.counts.and.with.oocyte.MII.FPKM.dt.gz")
{
    ##
    ## 5.1. summarize for each gene the spearman cor of RE-matched edit count (oocyte GV) ~ FPKM (oocyte MII)
    {
        ## compute spearman cor and count of patients
        GSE95477.patient.info.no.outliers.with.oocyte.GV.RE.matched.edit.counts.and.with.oocyte.MII.FPKM.dt -> .;
        .[, list(
            spearman.cor=cor(RE.matched.edits.count.oocyte.GV, FPKM.oocyte.MII, method="spearman"),
            count.of.patients=.N,
            RE.matched.edits.count.variable.type=c("Edit counts change", "Edit counts do not change")[(var(RE.matched.edits.count.oocyte.GV) == 0)+1]
        ), list(Gene_ID, Gene_Name, cluster)] -> .;
        ## select those observed in at least 6 patients to get a more reliable spearman.cor
        .[count.of.patients >= 6] -> .;
        ## prettify cluster
        .[, cluster.renamed.ordered:=factor(c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster])] -> .;
        ## RETURN
        .
    } -> ..to.plot.dt
    ##
    ## 5.2. plot distribution of spearman correlation coefficients (histogram)
    {
        ggplot(..to.plot.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=spearman.cor)) -> .;
        . + geom_histogram() -> .;
        . + theme_pubr() -> .;
        . + labs(x="Spearman correlation coefficient", y="# Genes", fill="") ->.;       
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.histogram.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##
    ## 5.3. plot distribution of spearman correlation coefficients for each cluster separately (density plot)
    {
        ggplot(..to.plot.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=spearman.cor, fill=cluster.renamed.ordered)) -> .;
        . + geom_density(alpha=0.5) -> .;
        . + theme_pubr() -> .;
        . + labs(x="Spearman correlation coefficient", y="density", fill="") ->.;
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.for.each.cluster.density.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##
    ## 5.4. plot distribution of spearman correlation coefficients for each cluster separately (boxplot)
    {
        ggplot(..to.plot.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=cluster.renamed.ordered, y=spearman.cor, fill=cluster.renamed.ordered)) -> .;
        . + geom_boxplot() -> .;
        . + stat_compare_means(comparisons=list(c("decay at 8-cell", "others")), method="wilcox.test", method.args=list(alternative="less")) ->.;
        . + theme_pubr() -> .;
        . + labs(x="", y="Spearman correlation coefficient", fill="") ->.;
        . + scale_y_continuous(limits=c(-1.0, 1.2)) -> .;
        . + theme(axis.text.x=element_blank()) -> .;
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.for.each.cluster.boxplot.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##
}

{
    ##
    ## 6.1. summarize for each gene the spearman cor of RE-matched edit count (oocyte GV) ~ FPKM (oocyte MII), in young and old separately
    {
        ## compute spearman cor and count of patients
        GSE95477.patient.info.no.outliers.with.oocyte.GV.RE.matched.edit.counts.and.with.oocyte.MII.FPKM.dt -> .;
        .[, list(
            spearman.cor=cor(RE.matched.edits.count.oocyte.GV, FPKM.oocyte.MII, method="spearman"),
            count.of.patients=.N,
            RE.matched.edits.count.variable.type=c("Edit counts change", "Edit counts do not change")[(var(RE.matched.edits.count.oocyte.GV) == 0)+1]
        ), list(Gene_ID, Gene_Name, cluster, maternal.bin.ordered)] -> .;
        ## select those observed in at least 3 patients to get a more reliable spearman.cor
        .[count.of.patients >= 3] -> .;
        ## prettify cluster
        .[, cluster.renamed.ordered:=factor(c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster])] -> .;
        ## RETURN
        .
    } -> ..to.plot.YAO.dt
    ##
    ## 6.2. plot distribution of spearman correlation coefficients (histogram)
    {
        ggplot(..to.plot.YAO.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=spearman.cor)) -> .;
        . + geom_histogram() -> .;
        . + theme_pubr() -> .;
        . + facet_grid(~maternal.bin.ordered) -> .;
        . + labs(x="Spearman correlation coefficient", y="# Genes", fill="") ->.;       
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.histogram.YAO.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##
    ## 6.3. plot distribution of spearman correlation coefficients for each cluster separately (density plot)
    {
        ggplot(..to.plot.YAO.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=spearman.cor, fill=maternal.bin.ordered)) -> .;
        . + geom_density(alpha=0.5) -> .;
        . + theme_pubr() -> .;
        . + facet_grid(cluster.renamed.ordered~.) -> .;
        . + labs(x="Spearman correlation coefficient", y="density", fill="") ->.;
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.for.each.cluster.density.YAO.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##
    ## 6.4. plot distribution of spearman correlation coefficients for each cluster separately (boxplot)
    {
        ggplot(..to.plot.YAO.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=maternal.bin.ordered, y=spearman.cor, fill=maternal.bin.ordered)) -> .;
        . + geom_boxplot() -> .;
        . + stat_compare_means(comparisons=list(c("old", "young")), method="wilcox.test", method.args=list(alternative="less")) ->.;
        . + theme_pubr() -> .;
        . + theme(axis.text.x=element_blank()) -> .;
        . + facet_grid(~cluster.renamed.ordered) -> .;
        . + labs(x="", y="Spearman correlation coefficient", fill="") ->.;
        . + scale_y_continuous(limits=c(-1.0, 1.2)) -> .;
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.for.each.cluster.boxplot.YAO.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##   
}




## Step 7. Examine for each gene RE-matched edit count ~ expression within the same sample
## NOTE: for each oocyte GV samples, only those genes with >=1 RE-matched edit were considered
## START

## RETURN
{
    ##
    ## 7.1. merge RE-matching edits count and expression for the same gene in the same sample
    ##
    {
        merge(
            x=GSE95477.RE.matched.edits.count.on.maternal.gene.only.dt[, list(SAMPLE, stage, maternal.bin.ordered, Gene_ID, Gene_Name, cluster, RE.matched.edits.count)],
            y=combined.gexpr.FPKM.melt.GSE95477.samples.only.dt[, list(SAMPLE=SAMPLE, Gene_ID, FPKM=FPKM)],
            by=c("Gene_ID", "SAMPLE"),
            all.x=FALSE, all.y=FALSE) -> .;
    } -> ..GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.dt;
    ##
    ## 7.2. filter against outlier samples
    ##
    {
        ..GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.dt[(SAMPLE %in% c("GSM2514781","GSM2514773")) == FALSE] -> .;
        .
    } -> ..GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.no.outlier.samples.dt;
    ##
    ## 7.2.BACKUP. write current result to file
    {
        fwrite(..GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.no.outlier.samples.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.no.outlier.samples.dt.gz")
    }
    ##
    ..GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.no.outlier.samples.dt
} -> GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.no.outlier.samples.dt



## Step 8. analyze the patterns 
## GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.no.outlier.samples.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.no.outlier.samples.dt.gz")
{
    ##
    ## 8.1. summarize for each gene the spearman cor of RE-matched edit count ~ FPKM within the same sample
    {
        ## compute spearman cor and count of patients
        GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.no.outlier.samples.dt -> .;
        .[, list(
            spearman.cor=cor(RE.matched.edits.count, FPKM, method="spearman"),
            count.of.patients=.N,
            RE.matched.edits.count.variable.type=c("Edit counts change", "Edit counts do not change")[(var(RE.matched.edits.count) == 0)+1]
        ), list(Gene_ID, Gene_Name, cluster, stage)] -> .;
        ## select those observed in at least 6 patients to get a more reliable spearman.cor
        .[count.of.patients >= 6] -> .;
        ## prettify cluster
        .[, cluster.renamed.ordered:=factor(c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster])] -> .;
        ## RETURN
        .
    } -> ..to.plot.within.sample.dt
    ##
    ## 8.2. plot distribution of spearman correlation coefficients (histogram) within same sample
    {
        ggplot(..to.plot.within.sample.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=spearman.cor)) -> .;
        . + geom_histogram() -> .;
        . + theme_pubr() -> .;
        . + facet_grid(~stage) ->.;
        . + labs(x="Spearman correlation coefficient", y="# Genes", fill="") ->.;       
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.histogram.within.sample.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##
    ## 8.3. plot distribution of spearman correlation coefficients for each cluster separately (density plot)
    {
        ggplot(..to.plot.within.sample.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=spearman.cor, fill=cluster.renamed.ordered)) -> .;
        . + geom_density(alpha=0.5) -> .;
        . + theme_pubr() -> .;
        . + facet_grid(~stage) ->.;
        . + labs(x="Spearman correlation coefficient", y="density", fill="") ->.;
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.for.each.cluster.density.within.sample.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##
    ## 8.4. plot distribution of spearman correlation coefficients for each cluster separately (boxplot)
    {
        ggplot(..to.plot.within.sample.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=cluster.renamed.ordered, y=spearman.cor, fill=cluster.renamed.ordered)) -> .;
        . + geom_boxplot() -> .;
        . + facet_grid(~stage) ->.;
        . + stat_compare_means(comparisons=list(c("decay at 8-cell", "others")), method="wilcox.test", method.args=list(alternative="less")) ->.;
        . + theme_pubr() -> .;
        . + labs(x="", y="Spearman correlation coefficient", fill="") ->.;
        . + scale_y_continuous(limits=c(-1.0, 1.2)) -> .;
        . + theme(axis.text.x=element_blank()) -> .;
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.for.each.cluster.boxplot.within.sample.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##
}


## Step 9. analyze the patterns 

{
    ##
    ## 9.1. summarize for each gene the spearman cor of RE-matched edit count ~ FPKM within the same sample, in (young and old, stage) separately
    {
        ## compute spearman cor and count of patients
        GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.no.outlier.samples.dt -> .;
        .[, list(
            spearman.cor=cor(RE.matched.edits.count, FPKM, method="spearman"),
            count.of.patients=.N,
            RE.matched.edits.count.variable.type=c("Edit counts change", "Edit counts do not change")[(var(RE.matched.edits.count) == 0)+1]
        ), list(Gene_ID, Gene_Name, cluster, maternal.bin.ordered, stage)] -> .;
        ## select those observed in at least 3 patients to get a more reliable spearman.cor
        .[count.of.patients >= 3] -> .;
        ## prettify cluster
        .[, cluster.renamed.ordered:=factor(c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster])] -> .;
        ## RETURN
        .
    } -> ..to.plot.within.sample.YAO.dt
    ##
    ## 9.2. plot distribution of spearman correlation coefficients (histogram)
    {
        ggplot(..to.plot.within.sample.YAO.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=spearman.cor)) -> .;
        . + geom_histogram() -> .;
        . + theme_pubr() -> .;
        . + facet_grid(~stage+maternal.bin.ordered) -> .;
        . + labs(x="Spearman correlation coefficient", y="# Genes", fill="") ->.;       
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.histogram.within.sample.YAO.png"), plot=., width.r=0.9, height.r=0.3)
    }
    ##
    ## 9.3. plot distribution of spearman correlation coefficients for each cluster separately (density plot)
    {
        ggplot(..to.plot.within.sample.YAO.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=spearman.cor, fill=maternal.bin.ordered)) -> .;
        . + geom_density(alpha=0.5) -> .;
        . + theme_pubr() -> .;
        . + facet_grid(cluster.renamed.ordered~stage) -> .;
        . + labs(x="Spearman correlation coefficient", y="density", fill="") ->.;
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.for.each.cluster.density.within.sample.YAO.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##
    ## 9.4. plot distribution of spearman correlation coefficients for each cluster separately (boxplot)
    {
        ggplot(..to.plot.within.sample.YAO.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=maternal.bin.ordered, y=spearman.cor, fill=maternal.bin.ordered)) -> .;
        . + geom_boxplot() -> .;
        . + stat_compare_means(comparisons=list(c("old", "young")), method="wilcox.test", method.args=list(alternative="greater")) ->.;
        . + theme_pubr() -> .;
        . + theme(axis.text.x=element_blank()) -> .;
        . + facet_grid(cluster.renamed.ordered~stage) -> .;
        . + labs(x="", y="Spearman correlation coefficient", fill="") ->.;
        . + scale_y_continuous(limits=c(-1.0, 1.4)) -> .;
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.for.each.cluster.boxplot.within.sample.YAO.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##
}



## 10. within-sample comparison, FPKM-thresholded
{
    ##
    ## 10.1. summarize for each gene the spearman cor of RE-matched edit count ~ FPKM within the same sample, in (young and old, stage) separately, only considering those with FPKM >= 5 in all samples
    {
        ## compute spearman cor and count of patients
        copy(GSE95477.RE.matched.edits.count.and.FPKM.on.maternal.gene.only.no.outlier.samples.dt) -> .;
        ## keep only those genes with FPKM >= 5 in all samples of the same stage
        .[, min.FPKM.in.the.current.stage:=min(FPKM), list(Gene_ID, stage)] -> .;
        .[min.FPKM.in.the.current.stage>=5] -> .;
        .[, list(
            spearman.cor=cor(RE.matched.edits.count, FPKM, method="spearman"),
            count.of.patients=.N,
            RE.matched.edits.count.variable.type=c("Edit counts change", "Edit counts do not change")[(var(RE.matched.edits.count) == 0)+1]
        ), list(Gene_ID, Gene_Name, cluster, maternal.bin.ordered, stage)] -> .;
        ## select those observed in at least 3 patients to get a more reliable spearman.cor
        .[count.of.patients >= 3] -> .;
        ## prettify cluster
        .[, cluster.renamed.ordered:=factor(c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster])] -> .;
        ## prettify stage
        temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
        merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] -> .;
        ## RETURN
        .
    } -> ..to.plot.within.sample.YAO.FPKM.thresholded.dt
    ##
    ## 10.2. plot distribution of spearman correlation coefficients for each cluster separately (boxplot)
    {
        ggplot(..to.plot.within.sample.YAO.FPKM.thresholded.dt[RE.matched.edits.count.variable.type == "Edit counts change"], aes(x=maternal.bin.ordered, y=spearman.cor, fill=maternal.bin.ordered)) -> .;
        . + geom_hline(yintercept=0, linetype="dashed", color="grey") -> .;
        . + geom_boxplot() -> .;
        . + stat_compare_means(comparisons=list(c("old", "young")), method="wilcox.test", method.args=list(alternative="greater")) ->.;
        . + theme_pubr() -> .;
        . + theme(axis.text.x=element_blank()) -> .;
        . + facet_grid(cluster.renamed.ordered~stage.description.ordered) -> .;
        . + labs(x="", y="Spearman's correlation coefficient", fill="") ->.;
        . + scale_y_continuous(limits=c(-1.0, 1.4)) -> .;
        ggsave.A4(filename=glue("{output.directory}/spearman.cor.for.each.cluster.boxplot.within.sample.YAO.FPKM.thresholded.png"), plot=., width.r=0.45, height.r=0.3)
    }
    ##
}


