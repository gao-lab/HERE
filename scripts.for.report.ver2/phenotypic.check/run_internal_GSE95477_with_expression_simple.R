library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("scales")
library("ggdendro")
library("GEOquery")
library("glue")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.RE.matching.edits.and.expression.simple/"
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
    ## ADD maternal bin
    .[, patient.maternal.bin:=c("AMA"="old", "YNG"="young")[sub(pattern=".*(YNG|AMA).*", replacement="\\1", x=title)]] -> .;
    ## ADD stage
    .[, stage:=c("developmental stage: germinal vesicle"="oocyte.GV", "developmental stage: metaphase II"="oocyte.MII")[developmental.stage]]
    ## KEEP COLUMNS SAMPLE, patient.ID
    .[, list(SAMPLE, patient.ID, patient.maternal.bin, stage)] -> .;
    ## RETURN
    .
} -> GSE95477.patient.melt.info.dt

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
         GSE95477.patient.melt.info.dt[, SAMPLE] -> .;
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
        fwrite(..combined.gexpr.FPKM.melt.GSE95477.samples.only.dt, glue("{output.directory}/combined.gexpr.FPKM.melt.GSE95477.samples.only.dt.gz"))
    }
    ##
    ## RETURN
    ..combined.gexpr.FPKM.melt.GSE95477.samples.only.dt
} -> combined.gexpr.FPKM.melt.GSE95477.samples.only.dt



## Step 3. Get RE-matched edit count __per MBS type__ in GSE95477 non-outlier samples
## START
## |<- 3.1. Get RE and relevant info present in GSE95477 stages [..GSE95477.RE.and.gene.per.stage.dt]
## 3.2. Filter all GSE95477 edits for those matching the REs above [..GSE95477.RE.matched.edits.dt]
## 3.3. Annotate MBS type per edit [..GSE95477.RE.matched.edits.with.MBS.annotation.dt]
## 3.4. Compute RE-matched edit count per (Gene, SAMPLE, MBS type) [..GSE95477.RE.matched.edits.with.MBS.annotation.count.dt]
## |<- 3.5. Get maternal gene info [..maternal.gene.info.dt]
## 3.6. Keep maternal genes only [..GSE95477.RE.matched.edits.with.MBS.annotation.count.on.maternal.gene.only.dt]
## 3.6.BACKUP.
##
{
    
    ## 3.1. Get RE and relevant info (i.e., position and their genes, and stage) present in GSE95477 stages (i.e.,oocytes GV and oocytes MII)
    {
        fread(
            "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz",
            select=c("CHROM", "POS", "Gene_ID", "Gene_Name", "stage")
        )[stage %in% c('oocyte.GV', 'oocyte.MII')] -> .;
        unique(.) -> .
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
    ## 3.3. Annotate MBS type per edit
    {
        fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/edited.ts.human.compared.with.original.annotated.summary.gene.and.edit.level.dt.csv.gz") -> .;
        merge(x=..GSE95477.RE.matched.edits.dt, y=.,
              by.x=c("POS", "Gene_ID"), by.y=c("edit.POS", "gene.id"),
              all.x=TRUE, all.y=FALSE) ->.;
    } -> ..GSE95477.RE.matched.edits.with.MBS.annotation.dt;
    ##
    ## 3.4. Compute RE-matched edit count per and MBS type (Gene, SAMPLE, MBS type)
    ## NOTE: about 1/3 genes do not have RE-matching edits that affect MBS
    {
        ..GSE95477.RE.matched.edits.with.MBS.annotation.dt[
           ,
            list(RE.matched.edits.count=.N),
            list(
                Gene_ID, Gene_Name,
                SAMPLE, gse, stage, maternal.age, maternal.bin.ordered=factor(c("young", "old")[(as.integer(maternal.age) > 35)+1], levels=c("old", "young")),
                gains.miRNA.sites)] -> .;
    } -> ..GSE95477.RE.matched.edits.with.MBS.annotation.count.dt;
    ##
    ## 3.5. Get maternal gene info
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
    ## 3.6. Keep maternal genes only
    {
        merge(x=..GSE95477.RE.matched.edits.with.MBS.annotation.count.dt, y=..maternal.gene.info.dt,
          by=c("Gene_ID", "Gene_Name"),
          all.x=FALSE, all.y=FALSE) ->.;
    } -> ..GSE95477.RE.matched.edits.with.MBS.annotation.count.on.maternal.gene.only.dt;
    ##
    
    ## 3.6.BACKUP. write current result to file
    {
        fwrite(..GSE95477.RE.matched.edits.with.MBS.annotation.count.on.maternal.gene.only.dt, glue("{output.directory}/GSE95477.RE.matched.edits.with.MBS.annotation.count.on.maternal.gene.only.dt.gz"))
    }
    
    ##
    ..GSE95477.RE.matched.edits.with.MBS.annotation.count.on.maternal.gene.only.dt
} -> GSE95477.RE.matched.edits.with.MBS.annotation.count.on.maternal.gene.only.dt; 

## Step 4.  filter for gene x stage  x MBS with drastic edit count difference
{

    ## 4.1. remove outlier samples
    {
        GSE95477.RE.matched.edits.with.MBS.annotation.count.on.maternal.gene.only.dt[(SAMPLE %in% c("GSM2514781","GSM2514773")) == FALSE] -> .;
    } -> ..GSE95477.RE.matched.edits.with.MBS.annotation.count.on.maternal.gene.only.no.outlier.samples.dt

    ## 4.2. compute info
    {
        ..GSE95477.RE.matched.edits.with.MBS.annotation.count.on.maternal.gene.only.no.outlier.samples.dt[, list(
            maximal.RE.matched.edits.count.in.old.patients=.SD[maternal.bin.ordered=='old', as.numeric(max(RE.matched.edits.count))],
            minimal.RE.matched.edits.count.in.young.patients=.SD[maternal.bin.ordered=='young', as.numeric(min(RE.matched.edits.count))],
            median.RE.matched.edits.count.in.old.patients=.SD[maternal.bin.ordered=='old', as.numeric(median(RE.matched.edits.count))],
            median.RE.matched.edits.count.in.young.patients=.SD[maternal.bin.ordered=='young', as.numeric(median(RE.matched.edits.count))],
            count.of.old.patients=sum(maternal.bin.ordered=='old'),
            count.of.young.patients=sum(maternal.bin.ordered=='young')
        ), list(Gene_ID, Gene_Name, cluster, stage, gains.miRNA.sites)] -> .;
    } -> ..GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.dt;

    ## 4.3. filter for genes with enough observations
    {
        ..GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.dt[(
            (stage == 'oocyte.GV' & count.of.old.patients >= 3 & count.of.young.patients >= 3) |
            (stage == 'oocyte.MII' & count.of.old.patients >= 4 & count.of.young.patients >= 4)
        )] -> .
    } -> ..GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.dt;

    ## 4.4. filter for genes with drastic RE-matched edits count difference
    {
        ..GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.dt[
            median.RE.matched.edits.count.in.young.patients - maximal.RE.matched.edits.count.in.old.patients > 0
        ] -> .;
    } -> ..GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.with.drastic.count.difference.dt
    
    ## 4.4.BACKUP
    {
        fwrite(..GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.with.drastic.count.difference.dt, glue("{output.directory}/GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.with.drastic.count.difference.dt.gz"))
    }

    ## 4.5. filter for genes with drastic RE-matched edits count difference (more strict)
    {
        ..GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.dt[
            minimal.RE.matched.edits.count.in.young.patients - maximal.RE.matched.edits.count.in.old.patients > 0
        ] -> .;
    } -> ..GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.with.drastic.count.difference.more.strict.dt
    
    ## 4.5.BACKUP
    {
        fwrite(..GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.with.drastic.count.difference.more.strict.dt, glue("{output.directory}/GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.with.drastic.count.difference.dt.more.strict.gz"))
    }

    ..GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.with.drastic.count.difference.dt
} -> GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.with.drastic.count.difference.dt;


## Step 5. compute expression change for each gene x stage
{
    ## 5.1. merge expression and sample info
    {
        merge(x=combined.gexpr.FPKM.melt.GSE95477.samples.only.dt, y=GSE95477.patient.melt.info.dt,
              by.x="SAMPLE", by.y="SAMPLE",
              all.x=TRUE, all.y=FALSE) -> .;
    } -> ..combined.gexpr.FPKM.melt.GSE95477.samples.only.with.sample.info.dt

    ## 5.2. discard outliers
    {
        ..combined.gexpr.FPKM.melt.GSE95477.samples.only.with.sample.info.dt[(SAMPLE %in%  c("GSM2514781","GSM2514773")) == FALSE] -> .;
    } -> ..combined.gexpr.FPKM.melt.GSE95477.samples.only.with.sample.info.no.outlier.samples.dt

    ## 5.2.BACKUP
    {
        fwrite(..combined.gexpr.FPKM.melt.GSE95477.samples.only.with.sample.info.no.outlier.samples.dt, glue("{output.directory}/combined.gexpr.FPKM.melt.GSE95477.samples.only.with.sample.info.no.outlier.samples.dt.gz"))
    }
    
    ..combined.gexpr.FPKM.melt.GSE95477.samples.only.with.sample.info.no.outlier.samples.dt
} ->  combined.gexpr.FPKM.melt.GSE95477.samples.only.with.sample.info.no.outlier.samples.dt


## Step 6. confirm whether the expression also falls
{
    
    ## 6.1. merge genes-of-interest and expression table
    {
        merge(x=combined.gexpr.FPKM.melt.GSE95477.samples.only.with.sample.info.no.outlier.samples.dt,
              y=GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.with.drastic.count.difference.dt,
              by.x=c("Gene_ID", "stage"), by.y=c("Gene_ID", "stage"),
              all.x=FALSE, all.y=TRUE) -> .;
    } -> ..gexpr.table.for.genes.whose.editing.drastically.dropped.in.old.GV.only.dt;

    ## 6.1.BACKUP.
    fwrite(..gexpr.table.for.genes.whose.editing.drastically.dropped.in.old.GV.only.dt, glue("{output.directory}/gexpr.table.for.genes.whose.editing.drastically.dropped.in.old.GV.only.dt.gz"))

    ## 6.2. plot
    {
        copy(..gexpr.table.for.genes.whose.editing.drastically.dropped.in.old.GV.only.dt) -> .;
        .[, list(FPKM.diff.log2.median.old.to.young=log2(.SD[patient.maternal.bin=='old', median(FPKM)]) - log2(.SD[patient.maternal.bin=='young', median(FPKM)])), list(Gene_ID, Gene_Name, cluster, stage, gains.miRNA.sites)] -> .;
        .[, gains.miRNA.sites.prettified:=c("does not gain MBS", "gains MBS")[gains.miRNA.sites + 1 ]][is.na(gains.miRNA.sites) == TRUE, gains.miRNA.sites.prettified:="no MBS identified\nupon editing"] -> .;
        . -> ..to.plot.dt;
        ggplot(..to.plot.dt, aes(x="a", y=FPKM.diff.log2.median.old.to.young)) -> .;
        . + geom_boxplot() -> .;
        ## . + stat_compare_means(comparisons=list(c("old", "young")), method="wilcox.test", method.args=list(alternative="less")) ->.;
        . + facet_grid(gains.miRNA.sites.prettified~stage+cluster ) -> .;
        ##. + labs(x="Spearman correlation coefficient", y="# Genes", fill="") ->.;       
        ggsave.A4(filename=glue("{output.directory}/FPKM.boxplot.between.young.and.old.on.editing.drastically.dropped.genes.png"), plot=., width.r=0.75, height.r=0.5)
    }
    
}



## Step 7. confirm whether the expression also falls
{
    
    ## 7.1. merge genes-of-interest and expression table
    {
        merge(x=combined.gexpr.FPKM.melt.GSE95477.samples.only.with.sample.info.no.outlier.samples.dt,
              y=GSE95477.maternal.bin.edit.count.info.per.gene.x.stage.x.MBS.annotation.type.with.enough.observations.with.drastic.count.difference.more.strict.dt,
              by.x=c("Gene_ID", "stage"), by.y=c("Gene_ID", "stage"),
              all.x=FALSE, all.y=TRUE) -> .;
    } -> ..gexpr.table.for.genes.whose.editing.drastically.dropped.in.old.GV.only.more.strict.dt;

    ## 7.1.BACKUP.
    fwrite(..gexpr.table.for.genes.whose.editing.drastically.dropped.in.old.GV.only.more.strict.dt, glue("{output.directory}/gexpr.table.for.genes.whose.editing.drastically.dropped.in.old.GV.only.more.strict.dt.gz"))

    ## 7.2. plot
    {
        copy(..gexpr.table.for.genes.whose.editing.drastically.dropped.in.old.GV.only.more.strict.dt) -> .;
        .[, list(FPKM.diff.log2.median.old.to.young=log2(.SD[patient.maternal.bin=='old', median(FPKM)]) - log2(.SD[patient.maternal.bin=='young', median(FPKM)])), list(Gene_ID, Gene_Name, cluster, stage, gains.miRNA.sites)] -> .;
        .[, gains.miRNA.sites.prettified:=c("does not gain MBS", "gains MBS")[gains.miRNA.sites + 1 ]][is.na(gains.miRNA.sites) == TRUE, gains.miRNA.sites.prettified:="no MBS identified\nupon editing"] -> .;
        . -> ..to.plot.dt;
        ggplot(..to.plot.dt, aes(x="a", y=FPKM.diff.log2.median.old.to.young)) -> .;
        . + geom_boxplot() -> .;
        ## . + stat_compare_means(comparisons=list(c("old", "young")), method="wilcox.test", method.args=list(alternative="less")) ->.;
        . + facet_grid(gains.miRNA.sites.prettified~stage+cluster ) -> .;
        ##. + labs(x="Spearman correlation coefficient", y="# Genes", fill="") ->.;       
        ggsave.A4(filename=glue("{output.directory}/FPKM.boxplot.between.young.and.old.on.editing.drastically.dropped.genes.more.strict.png"), plot=., width.r=0.75, height.r=0.5)
    }
    
}



