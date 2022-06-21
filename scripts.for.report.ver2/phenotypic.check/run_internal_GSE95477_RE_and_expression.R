library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("scales")
library("GEOquery")
source("./scripts/common/ggpubr.A4.R")

## read total edits in all samples
edit.info.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")

## read other tables (for filtering)
## 1. RE and their genes
subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")
## 2. annotation of maternal genes
combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt <- fread("result/S42_1__annotate_embryonic_genes/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/201221-fifth-phenotype-collection/combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt.txt.gz")
## 3. original expression table
combined.gexpr.FPKM.dt <- {
    read.table("result/BS06_1__get_expression_level/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/combined.gexpr.FPKM.matrix.txt") -> .;
    data.table(Gene_ID=rownames(.), .) -> .;
    .
}

combined.gexpr.FPKM.melt.dt <- melt(combined.gexpr.FPKM.dt, id.vars="Gene_ID", variable.name="SAMPLE", value.name="FPKM")

normal.RE.and.gene.per.stage.dt <- {
    subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt ->.;
    ##s keep matched stages only
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell')] ->.;
    ##= keep columns needed 
    .[, list(CHROM, POS, Gene_ID, Gene_Name, Annotation.corrected, stage)] %>% unique ->.;
}

maternal.gene.info.dt <- combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt[cluster %in% c("maternal.decay", "maternal.others"), list(Gene_ID=gene.id, Gene_Name=gene.name, cluster)]

edit.normal.RE.only.maternal.gene.only.dt <- {
    edit.info.dt ->.;
    ##s keep matched stages only
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', 'zygote.2PN', '2-cell', '4-cell', '8-cell')] ->.;
    .[, stage.corrected:=stage][stage=="zygote.2PN", stage.corrected:="zygote"] ->.;
    ##m keep normal-RE edits only
    merge(x=., y=normal.RE.and.gene.per.stage.dt,
          by.x=c("CHROM", "POS", "stage.corrected"), by.y=c("CHROM", "POS", "stage"),
          all.x=FALSE, all.y=FALSE, allow.cartesian=TRUE) ->.;
    ## keep maternal genes only
    merge(x=., y=maternal.gene.info.dt,
          by=c("Gene_ID", "Gene_Name"),
          all.x=FALSE, all.y=FALSE) ->.;
}

GSE95477.spearman.cor.between.normal.RE.count.and.FPKM.across.Gene_IDs.per.nonoutlier.SAMPLE.dt <- {
    edit.normal.RE.only.maternal.gene.only.dt -> .;
    ## select GSE95477 samples only
    .[gse=="GSE95477"] -> .;
    ## FILTER ROWS against outlier SAMPLE's: SAMPLE %in% c("GSM2514781", "GSM2514773") == FALSE
    .[SAMPLE %in% c("GSM2514781", "GSM2514773") == FALSE] ->.;
    ## SUMMARIZE ACROSS ID's FOR EACH ~SAMPLE(with gse, stage, maternal.age, maternal.bin.ordered) + Gene_ID(wtih Gene_Name, cluster): normal.RE.count=.N
    .[, list(normal.RE.count=.N), list(Gene_ID, Gene_Name, cluster, gse, SAMPLE, stage, maternal.age, maternal.bin.ordered=factor(c("young", "old")[(as.integer(maternal.age) > 35)+1], levels=c("old", "young")))] -> .;
    ## WRITE current reult (SAMPLE, Gene_ID, normal.RE.count) to file
    {
        fwrite(., "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.normal.RE.count.per.Gene_IDs.and.nonoutlier.SAMPLE.dt.gz")
    }
    ## LEFT-JOIN with combined.gexpr.FPKM.melt.dt ON ~Gene_ID + SAMPLE
    merge(x=., y=combined.gexpr.FPKM.melt.dt, by=c("Gene_ID", "SAMPLE"), all.x=TRUE, all.y=FALSE) -> .;
    ## WRITE current result to file
    {
        fwrite(., "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/temp.GSE95477.normal.RE.count.and.FPKM.across.Gene_IDs.per.nonoutlier.SAMPLE.dt.txt.gz")
    }
    ## SUMMARIZE ACROSS GENE_ID's FOR EACH ~SAMPLE(with stage, cluster, maternal.bin.ordered): spearman.cor=cor(normal.RE.count, FPKM, method="spearman")
    .[, list(spearman.cor=cor(normal.RE.count, FPKM, method="spearman")), list(SAMPLE, stage, cluster, maternal.bin.ordered)] -> .;
    ## RETURN
    .
}

## fwrite(GSE95477.spearman.cor.between.normal.RE.count.and.FPKM.across.Gene_IDs.per.nonoutlier.SAMPLE.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.spearman.cor.between.normal.RE.count.and.FPKM.across.Gene_IDs.per.nonoutlier.SAMPLE.dt.txt.gz")

## diagnosis 1
## for each gene, oocyte.GV normal RE count ~ oocyte.MII (same patient) FPKM

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
    ## RETURN
    .
} -> GSE95477.patient.info.dt

## Step 2. Get expression tables for GSE95477
## NOTE: do not use the edit x FPKM table, because it does not contain genes with no edits per sample
{
    ## GET GSM accessions of GSE95477
    c(GSE95477.patient.info.dt[, oocyte.GV], GSE95477.patient.info.dt[, oocyte.MII]) -> .;
    ## FILTER ROW of (Gene_ID, SAMPLE, FPKM) for SAMPLE %in% GSE95477 samples
    combined.gexpr.FPKM.melt.dt[SAMPLE %in% .] -> .;
    ## WRITE current result to file
    fwrite(., "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/temp.GSE95477.exp.melt.dt.gz")
    ## RETURN
    .
} -> GSE95477.exp.melt.dt
    
## Step 3. Get 

GSE95477.normal.RE.count.per.Gene_IDs.and.nonoutlier.SAMPLE.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.normal.RE.count.per.Gene_IDs.and.nonoutlier.SAMPLE.dt.gz")

{
    ## GET columns of GSE95477.patient.info.dt
    GSE95477.patient.info.dt[, list(patient.ID, SAMPLE.oocyte.GV=oocyte.GV, SAMPLE.oocyte.MII=oocyte.MII)] -> .;
    ## FILTER ROW against outlier samples
    .[(SAMPLE.oocyte.GV %in% c("GSM2514781", "GSM2514773") == FALSE) & (SAMPLE.oocyte.MII %in% c("GSM2514781", "GSM2514773") == FALSE)] -> .;
    ## LEFT-JOIN with (SAMPLE, Gene_ID, normal.RE.count) ON SAMPLE.oocyte.GV==SAMPLE
    merge(x=., y=GSE95477.normal.RE.count.per.Gene_IDs.and.nonoutlier.SAMPLE.dt[, list(SAMPLE.oocyte.GV=SAMPLE, maternal.bin.ordered, Gene_ID, Gene_Name, cluster, normal.RE.count.oocyte.GV=normal.RE.count)],
          by="SAMPLE.oocyte.GV",
          all.x=TRUE, all.y=FALSE) -> .;
    ## LEFT-JOIN with (SAMPLE, Gene_ID, FPKM) ON (SAMPLE.oocyte.MII==oocyte.MII, Gene_ID)
    merge(x=., y=GSE95477.exp.melt.dt[, list(SAMPLE.oocyte.MII=SAMPLE, Gene_ID, FPKM.oocyte.MII=FPKM)],
          by=c("SAMPLE.oocyte.MII", "Gene_ID"),
          all.x=TRUE, all.y=FALSE) -> .;
    ## WRITE current result (Gene_ID, SAMPLE.oocyte.GV, normal.RE.count in GV, SAMPLE.oocyte.MII, FPKM in MII ) to file
    {
        fwrite(., "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.oocyte.GV.normal.RE.count.and.oocyte.MII.FPKM.per.Gene_ID.per.patient.dt.gz")
    }
} -> GSE95477.oocyte.GV.normal.RE.count.and.oocyte.MII.FPKM.per.Gene_ID.per.patient.dt

## TODO analyze the patterns of report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.oocyte.GV.normal.RE.count.and.oocyte.MII.FPKM.per.Gene_ID.per.patient.dt.gz
## GSE95477.oocyte.GV.normal.RE.count.and.oocyte.MII.FPKM.per.Gene_ID.per.patient.dt <- fread("report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.oocyte.GV.normal.RE.count.and.oocyte.MII.FPKM.per.Gene_ID.per.patient.dt.gz")


{
    GSE95477.oocyte.GV.normal.RE.count.and.oocyte.MII.FPKM.per.Gene_ID.per.patient.dt -> .;
    ## PLOT
    {
        ## ggplot(., aes(x=normal.RE.count.oocyte.GV, y=FPKM.oocyte.MII)) -> .;
        ## . + geom_point() -> .;
        ## ggplot(., aes(x=cut(normal.RE.count.oocyte.GV, breaks=c(0, 1, 2, 5, 10, 20, 50, Inf)), y=FPKM.oocyte.MII)) -> .;
        ggplot(., aes(x=cut(FPKM.oocyte.MII, breaks=c(0, 0.1, 1, 10, 100, 1000, Inf)), y=normal.RE.count.oocyte.GV, color=maternal.bin.ordered)) -> .;
        . + geom_boxplot() -> .;
        ## . + facet_grid(maternal.bin.ordered ~ patient.ID ~ cluster) -> .;
        ## . + facet_grid(maternal.bin.ordered ~ cluster) -> .;
        . + facet_grid(~ cluster) -> .;
        . + theme_pubr() -> .;
        ## . + scale_y_continuous(limits=c(0, 40)) ->.;
    }
    .
}

## TODO 从图上看，似乎RE越少，表达量越高，但是本来表达量高的基因就更容易成为RNA编辑的目标（假设所有RNA编辑发生概率均等的前提下），所以这里看到的规律不能直接归咎于GV期RE对MII期表达量的影响（他可能收到GV期表达量的影响），不能直接说GV期RE变少会导致MII期表达量变低
## TODO，但是我们可以看到，只要是decay基因，图上无论是哪个FPKM窗口，young的RE count数总是高于old的一些，说明这个现象在各种oocyte MII高表达模式的基因上都存在，是个普遍现象
## TODO 另一个和表达量的关系的点是，RE个数差距最大的地方在表达量中下的地方，暗示他可能主要作用在需要母源降解，但是MII期表达又不是很高的基因上

## diagnosis 2
temp.GSE95477.normal.RE.count.and.FPKM.across.Gene_IDs.per.nonoutlier.SAMPLE.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/temp.GSE95477.normal.RE.count.and.FPKM.across.Gene_IDs.per.nonoutlier.SAMPLE.dt.txt.gz")

{
    temp.GSE95477.normal.RE.count.and.FPKM.across.Gene_IDs.per.nonoutlier.SAMPLE.dt -> .;
    ## FILTER ROW
    ggplot(.[stage=='oocyte.GV']) -> .;
    ##. + geom_point(aes(x=normal.RE.count, y=FPKM)) -> .;
    ##. + geom_histogram(aes(x=normal.RE.count)) -> .;
    . + geom_boxplot(aes(x=SAMPLE, y=normal.RE.count, color=cut(FPKM, c(0, 0.1, 1, 5, 10, 50, 100, 500, 1000, Inf)))) + scale_y_continuous(breaks=seq(0, 65, 5)) + facet_grid(cluster + maternal.bin.ordered~., scales="free_y") + coord_flip() -> .;
    ## . + facet_wrap(cluster + maternal.bin.ordered~SAMPLE, scales="fixed") -> .;
    . + theme_pubr() -> .;
    .
}


## compute per sample: count of REs per gene
RE.count.info.dt <- edit.normal.RE.only.maternal.gene.only.dt[, list(count.RE=length(unique(ID)), count.gene=length(unique(Gene_ID)), mean.RE.per.gene=length(unique(ID)) / length(unique(Gene_ID)) ), list(SAMPLE, gse, stage, stage.corrected, is.normal, treatment, disease, maternal.age, cluster)]


## Plot all samples of GSE95477
{
    RE.count.info.dt[gse=="GSE95477"] ->.;
    ## bin maternal ages
    .[, maternal.bin.ordered:=factor(c("young", "old")[(as.integer(maternal.age) > 35)+1], levels=c("old", "young")) ]
    ##= prettify cluster
    .[, cluster.renamed.ordered:=factor(c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster])]
    ## prettify stages
    temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## start plotting
    ggplot(., aes(x=maternal.bin.ordered, y=mean.RE.per.gene, fill=maternal.bin.ordered)) ->.;
    ## add mean comparison
    . + stat_compare_means(comparisons=list(c("old", "young")), method="wilcox.test", method.args=list(alternative="less")) ->.;
    ## plot boxplots
    . + geom_boxplot() ->.;
    ## facet
    . + facet_grid(cluster.renamed.ordered~stage.description.ordered, scales="free_y") ->.;
    ## add theme
    . + theme_pubr() ->.;
    . + theme(axis.text.x = element_blank()) ->.;
    . + scale_y_continuous(limits=c(4.5, 6.5)) ->.;
    . + labs(x="", y="Average # RE-matching edits per gene", fill="Maternal age group") ->.;
    ## save image
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/RE.count.boxplot.GSE95477.png",
        plot=.,
        width.r=0.45, height.r=0.4)
}


## Plot clustering result

{
    edit.normal.RE.only.maternal.gene.only.dt ->.;
    ##s select GSE95477 only
    .[gse %in% c("GSE95477")] ->.;
    ##c get RE (on maternal genes) per sample and cluster
    .[, list(SAMPLE, stage, maternal.age, maternal.bin.ordered=factor(c("young", "old")[(as.integer(maternal.age) > 35)+1], levels=c("old", "young")))] %>% unique -> temp.phenotype.dt
    .[, list(CHROM, POS, SAMPLE, one=1)] %>% unique ->.;
    ##d dcast
    dcast(., CHROM + POS ~ SAMPLE, value.var="one", fill=0) ->.;
    ##t transpose
    t(.[, c(-1, -2)]) -> .;
    ## generate cluster plot
    dendro_data(hclust(dist(., method="manhattan"))) ->.;
    .[["labels"]] <- merge(x=.[["labels"]], y=temp.phenotype.dt, by.x="label", by.y="SAMPLE", all=TRUE)
    ggplot() + geom_segment(data=segment(.), aes(x=x, y=y, xend=xend, yend=yend)) + geom_text(data=label(.), aes(x=x, y=y, label=paste(sep="", label, " ; age: ", maternal.age), hjust=0, color=paste(sep="", stage, ";", maternal.bin.ordered))) + coord_flip() +  scale_y_continuous(expand=expansion(mult=c(0.3, 1.3)), breaks=seq(0, 3000, 1000), trans=reverse_trans()) -> .;
    . + theme_pubr() ->.;
    . + theme(axis.text.y=element_blank(), axis.ticks.y = element_blank()) ->.;
    . + guides(color=guide_legend(nrow=2)) -> .;
    . + labs(x="", y="hclust height", color="") ->.;
    ggsave.A4(
        "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/GSE95477.oocyte.hclust.png",
        .,
        width.r=0.45, height.r=0.4
    )
}

## Plot all samples of GSE95477 except for "GSM2514781","GSM2514773"
{
    RE.count.info.dt[gse=="GSE95477"] ->.;
    ## remove GSM2514781
    .[SAMPLE %in% c("GSM2514781","GSM2514773")==FALSE] ->.;
    ## bin maternal ages
    .[, maternal.bin.ordered:=factor(c("young", "old")[(as.integer(maternal.age) > 35)+1], levels=c("old", "young")) ]
    ##= prettify cluster
    .[, cluster.renamed.ordered:=factor(c("maternal.decay"="decay at 8-cell", "maternal.others"="others")[cluster])]
    ## prettify stages
    temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## start plotting
    ggplot(., aes(x=maternal.bin.ordered, y=mean.RE.per.gene, fill=maternal.bin.ordered)) ->.;
    ## add mean comparison
    . + stat_compare_means(comparisons=list(c("old", "young")), method="wilcox.test", method.args=list(alternative="less")) ->.;
    ## plot boxplots
    . + geom_boxplot() ->.;
    ## facet
    . + facet_grid(cluster.renamed.ordered~stage.description.ordered, scales="free_y") ->.;
    ## add theme
    . + theme_pubr() ->.;
    . + theme(axis.text.x = element_blank()) ->.;
    . + scale_y_continuous(limits=c(4.5, 6.5)) ->.;
    . + labs(x="", y="Average # RE-matching edits per gene", fill="Maternal age group") ->.;
    ## save image
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/RE.count.boxplot.GSE95477.without.GSM2514781.and.GSM2514773.png",
        plot=.,
        width.r=0.45, height.r=0.3)
}


