library("data.table")
library("magrittr")
library("readxl")
library("foreach")
library("iterators")
library("scales")
source("./scripts/common/ggpubr.A4.R")

RE.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

edited.ts.human.compared.with.original.annotated.summary.gene.and.edit.level.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/edited.ts.human.compared.with.original.annotated.summary.gene.and.edit.level.dt.csv.gz")

MTC.info.dt <- fread("result/S42_1__annotate_embryonic_genes/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/201221-fifth-phenotype-collection/combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt.txt.gz")




temp.tables.list <- list()

{
    copy(RE.valid.genes.only.dt) ->.;
    ##s select 3'-UTR only
    .[grepl("3_prime_UTR_variant", Annotation.pasted)== TRUE] ->.;
    ##s select detected edits only per sample
    .[is.na(AC)==FALSE] ->.;
    ## merge with edit miRNA info
    merge(x=., y=edited.ts.human.compared.with.original.annotated.summary.gene.and.edit.level.dt,
          by.x=c("Gene_ID", "POS"), by.y=c("gene.id", "edit.POS"),
          all.x=FALSE, all.y=FALSE) ->.;
    ##c compute count of RE per edit-type x gene x sample
    .[, list(
        count.of.edits.gaining.miRNA.sites=sum(gains.miRNA.sites),
        count.of.edits.losing.miRNA.sites=sum(loses.miRNA.sites)
    ), list(Gene_ID, SAMPLE, stage)] ->.;
    ##save temp table
    temp.tables.list[["RE.count.per.gene.and.sample.and.stage.dt"]] <- copy(.)
    ## start plotting
    ##s select stage of interest only
    .[stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell")] -> .;
    ##melt data for plotting
    melt(data=., id.vars=c("Gene_ID", "SAMPLE", "stage"), measure.vars=c("count.of.edits.gaining.miRNA.sites", "count.of.edits.losing.miRNA.sites"), variable.name="metric", value.name="value") -> .;
    ##= prettify metrics
    .[, metric.renamed:=c("count.of.edits.gaining.miRNA.sites"="MBS gain", "count.of.edits.losing.miRNA.sites"="MBS lost")[metric]]
    ##= prettify values by clipping them
    .[, value.clipped:=pmin(value, 10)]
    ##= prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## copy a data for the statistical test
    temp.data.dt <- copy(.)
    ## run ggplot
    ggplot(., aes(x=metric.renamed, y=value.clipped)) -> .;
    . + geom_boxplot() ->.;
    . + stat_compare_means(data=temp.data.dt, mapping=aes(y=value), method="wilcox.test", paired=TRUE, comparisons=list(c("MBS gain", "MBS lost")), method.args=list(alternative="greater")) -> .;
    . + labs(x="effect of RE on MBS", y="# REs per (gene, sample)\n(upper-bounded by 10)") -> .;
    . + theme_pubr() ->.;
    . + scale_y_continuous(limits=c(0,12), breaks=seq(0, 10, 2)) ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/RE.MBS.gain.and.loss.count.boxplot.png",
        plot=.,
        width.r=0.4, height.r=0.25
    )
}


## plot per gene x sample the counts of MBS-gaining RE between MTC and other maternal genes
{
    copy(temp.tables.list[["RE.count.per.gene.and.sample.and.stage.dt"]]) ->.;
    ##m merge with cluster annotation
    merge(x=., y=MTC.info.dt[, list(Gene_ID=gene.id, cluster)],
          by="Gene_ID",
          all.x=TRUE, all.y=FALSE) ->.;
    ##s select maternal genes only
    .[cluster %in% c("maternal.decay", "maternal.others")] ->.;
    ##s select stages on interest only
    .[stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell")] -> .;
    ##= prettify clusters
    .[, cluster.renamed:=c("maternal.decay"="decay at\n8-cell", "maternal.others"="others")[cluster]]
    ##= prettify values by cutting them
    .[,
      count.of.edits.gaining.miRNA.sites.cut:=
          cut(count.of.edits.gaining.miRNA.sites, breaks=c(-Inf, 0:2, Inf), include.lowest=FALSE) %>%
          {temp.levels <- c("(-Inf,0]"="0", "(0,1]"="1", "(1,2]"="2", "(2,Inf]"=">=3"); factor(temp.levels[.], levels=temp.levels)}
      ]
    ##= prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## copy a data for the statistical test
    temp.data.dt <- copy(.)
    ## run statistical test
    temp.p.value.between.MTC.and.others <- wilcox.test(
        x=temp.data.dt[cluster=="maternal.decay", count.of.edits.gaining.miRNA.sites],
        y=temp.data.dt[cluster=="maternal.others", count.of.edits.gaining.miRNA.sites],
        alternative="greater"
    )$p.value %>% {ifelse(.<2.22e-16, "p < 2.22e-16", scientific(.))}
    temp.p.value.between.MTC.and.baseline <- wilcox.test(
        x=temp.data.dt[cluster=="maternal.decay", count.of.edits.gaining.miRNA.sites],
        y=NULL, mu=0,
        alternative="greater"
    )$p.value %>% {ifelse(.<2.22e-16, "p < 2.22e-16", scientific(.))}
    temp.stat.dt <- data.table(
        count.of.edits.gaining.miRNA.sites.cut="0",
        `.y.`="count",
        group1=rep("decay at\n8-cell", 2),
        group2=c("others", "baseline"),
        y.position=c(110, 125),
        p.value.raw=c(temp.p.value.between.MTC.and.others, temp.p.value.between.MTC.and.baseline)
    )
    ## create count summary for plot
    .[, list(count=.N), list(cluster.renamed, count.of.edits.gaining.miRNA.sites.cut)] -> .;
    .[, percentage:=count/sum(count) * 100, list(cluster.renamed)]
    rbindlist(list(., data.table(count=1, percentage=100, cluster.renamed="baseline", count.of.edits.gaining.miRNA.sites.cut="0")), use.names=TRUE) -> .;
    ## re-prettify cluster
    .[, cluster.renamed.ordered:=factor(cluster.renamed, c("decay at\n8-cell", "others", "baseline"))]
    ## run ggplot
    ggplot(., aes(x=cluster.renamed.ordered, y=percentage, fill=count.of.edits.gaining.miRNA.sites.cut)) -> .;
    . + geom_bar(stat="identity", position="stack") ->.;
    . + stat_pvalue_manual(data=temp.stat.dt, label="p.value.raw") -> .;
    . + labs(x="", y="% (gene, sample)", fill="# MBS-\ngaining\nREs on\nthe gene") -> .;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="right") ->.;
    . + scale_fill_manual(values=c("grey50", "#DD5532", "#FF6F39", "#FFA338")) ->.;
    . + scale_y_continuous(limits=c(0, 130), breaks=seq(0, 100, 25)) -> .;
    . + guides(fill=guide_legend(ncol=1)) ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/RE.MBS.gain.count.on.MTC.or.nonMTC.maternal.genes.barplot.png",
        plot=.,
        width.r=0.45, height.r=0.25
    )
}

##plot per gene x sample the counts of different MBS-net-change-type RE between MTC and other maternal genes
#### the filled bar version
{
    copy(temp.tables.list[["RE.count.per.gene.and.sample.and.stage.dt"]]) ->.;
    ##= add diff change info
    .[, `:=`(
        diff.between.count.of.edits.gaining.miRNA.sites.and.count.of.edits.losing.miRNA.sites=count.of.edits.gaining.miRNA.sites - count.of.edits.losing.miRNA.sites
    )] ->.;
    ##m merge with cluster annotation
    merge(x=., y=MTC.info.dt[, list(Gene_ID=gene.id, cluster)],
          by="Gene_ID",
          all.x=TRUE, all.y=FALSE) ->.;
    ##s select maternal genes only
    .[cluster %in% c("maternal.decay", "maternal.others")] ->.;
    ##s select stages on interest only
    .[stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell")] -> .;
    ##= prettify clusters
    .[, cluster.renamed:=c("maternal.decay"="decay at\n8-cell", "maternal.others"="others")[cluster]]
    ##= prettify values by cutting them
    .[,
      diff.between.count.of.edits.gaining.miRNA.sites.and.count.of.edits.losing.miRNA.sites.cut:=
          cut(diff.between.count.of.edits.gaining.miRNA.sites.and.count.of.edits.losing.miRNA.sites, breaks=c(-Inf, -3:2, Inf), include.lowest=FALSE) %>%
          {temp.levels <- c("(-Inf,-3]"="<=-3", "(-3,-2]"="-2", "(-2,-1]"="-1", "(-1,0]"="0", "(0,1]"="1", "(1,2]"="2", "(2,Inf]"=">=3"); factor(temp.levels[.], levels=temp.levels)}
      ]
    ##= prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## copy a data for the statistical test
    temp.data.dt <- copy(.)
    ## run statistical test
    temp.p.value.between.MTC.and.others <- wilcox.test(
        x=temp.data.dt[cluster=="maternal.decay", diff.between.count.of.edits.gaining.miRNA.sites.and.count.of.edits.losing.miRNA.sites],
        y=temp.data.dt[cluster=="maternal.others", diff.between.count.of.edits.gaining.miRNA.sites.and.count.of.edits.losing.miRNA.sites],
        alternative="greater"
    )$p.value %>% {ifelse(.<2.22e-16, "p < 2.22e-16", scientific(.))}
    temp.p.value.between.MTC.and.baseline <- wilcox.test(
        x=temp.data.dt[cluster=="maternal.decay", diff.between.count.of.edits.gaining.miRNA.sites.and.count.of.edits.losing.miRNA.sites],
        y=NULL, mu=0,
        alternative="greater"
    )$p.value %>% {ifelse(.<2.22e-16, "p < 2.22e-16", scientific(.))}
    temp.stat.dt <- data.table(
        diff.between.count.of.edits.gaining.miRNA.sites.and.count.of.edits.losing.miRNA.sites.cut="0",
        `.y.`="count",
        group1=rep("decay at\n8-cell", 2),
        group2=c("others", "baseline"),
        y.position=c(110, 125),
        p.value.raw=c(temp.p.value.between.MTC.and.others, temp.p.value.between.MTC.and.baseline)
    )
    ## create count summary for plot
    .[, list(count=.N), list(cluster.renamed, diff.between.count.of.edits.gaining.miRNA.sites.and.count.of.edits.losing.miRNA.sites.cut)] -> .;
    .[, percentage:=count/sum(count) * 100, list(cluster.renamed)]
    rbindlist(list(., data.table(count=1, percentage=100, cluster.renamed="baseline", diff.between.count.of.edits.gaining.miRNA.sites.and.count.of.edits.losing.miRNA.sites.cut="0")), use.names=TRUE) -> .;
    ## re-prettify cluster
    .[, cluster.renamed.ordered:=factor(cluster.renamed, c("decay at\n8-cell", "others", "baseline"))]
    ## run ggplot
    ggplot(., aes(x=cluster.renamed.ordered, y=percentage, fill=diff.between.count.of.edits.gaining.miRNA.sites.and.count.of.edits.losing.miRNA.sites.cut)) -> .;
    . + geom_bar(stat="identity", position="stack") ->.;
    . + stat_pvalue_manual(data=temp.stat.dt, label="p.value.raw") -> .;
    . + labs(x="", y="% (gene, sample)", fill="Net change of MBS\nby RE on the gene\n(defined as\n# MBS-gaining RE -\n# MBS-losing RE)") -> .;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top") ->.;
    . + scale_fill_manual(values=c("#69D2FF", "#00B6EB", "#0091FF", "grey50", "#DD5532", "#FF6F39", "#FFA338")) ->.;
    . + scale_y_continuous(limits=c(0, 130), breaks=seq(0, 100, 25)) -> .;
    . + guides(fill=guide_legend(nrow=4)) ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/RE.MBS.net.change.count.on.MTC.or.nonMTC.maternal.genes.barplot.png",
        plot=.,
        width.r=0.45, height.r=0.5
    )
}


