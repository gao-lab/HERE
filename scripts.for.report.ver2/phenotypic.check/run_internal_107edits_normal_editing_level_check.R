library("data.table")
library("readxl")
library("glue")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/lost.edits.coverage.check/"
dir.create(output.directory, recursive=TRUE)






## per-stage count, barplot
{
    

    ## 1. get the 107 edits
    {
        "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/disease.or.old.mother.embryo.lost.RE.dt.csv.gz" -> .;
        fread(.) -> .;
        unique(.[, list(CHROM, POS, stage, completely.lost.in.GSE95477.old.mother.embryos, completely.lost.in.GSE133854.AG.embryos, completely.lost.in.GSE133854.PG.embryos)]) -> .;
    } -> ..all.107edits.sites.dt

    ## 2. load their editing level in normal embryos

    {
        
        "result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz" -> .;
        fread(.) -> .;
        
        merge(
            x=.[is.normal==TRUE], y=..all.107edits.sites.dt,
            by=c("CHROM", "POS", "stage"),
            all.x=FALSE, all.y=TRUE
        ) -> .;

        .
    } -> ..all.107edits.normal.editing.level.dt

    ## 3. plot
    {
        
        copy(..all.107edits.normal.editing.level.dt) -> .;
        ## prettify stage
        temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
        merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
        ## prettify label
        .[completely.lost.in.GSE133854.AG.embryos == TRUE, label:="AG"] -> .;
        .[completely.lost.in.GSE133854.PG.embryos == TRUE, label:="PG"] -> .;
        .[completely.lost.in.GSE95477.old.mother.embryos == TRUE, label:="elder mother"] -> .;
        ##
        . -> ..to.plot.dt

        ## plot editing level (GSE95477 elder mother lost)
        ggplot(..to.plot.dt[completely.lost.in.GSE95477.old.mother.embryos == TRUE], aes(x=paste(sep="", CHROM, "_", POS), y=AF)) ->.;
        . + geom_boxplot() ->.;
        . + coord_flip() -> .;
        . + theme_pubr() ->.;
        . + theme(strip.text.y.right = element_text(angle = 0)) -> .;
        . + facet_grid(stage.description.ordered~., scales="free_y", space="free_y") -> .;
        . + labs(x="", y="Editing level") ->.;
        . + scale_y_continuous(limits=c(0.1, 1), breaks=seq(0.1, 1, 0.3)) -> .;
        ggsave.A4(
            filename=glue("{output.directory}/all.107edits.normal.editing.level.GSE95477.boxplot.png"),
            plot=.,
            width.r=0.9, height.r=0.9
        )

        ## plot editing level (GSE133854 AG lost)
        ggplot(..to.plot.dt[completely.lost.in.GSE133854.AG.embryos == TRUE], aes(x=paste(sep="", CHROM, "_", POS), y=AF)) ->.;
        . + geom_boxplot() ->.;
        . + coord_flip() -> .;
        . + theme_pubr() ->.;
        . + theme(strip.text.y.right = element_text(angle = 0)) -> .;
        . + facet_grid(stage.description.ordered ~., scales="free_y", space="free_y") -> .;
        . + labs(x="", y="Editing level") ->.;
        . + scale_y_continuous(limits=c(0.1, 1), breaks=seq(0.1, 1, 0.3)) -> .;
        ggsave.A4(
            filename=glue("{output.directory}/all.107edits.normal.editing.level.GSE133854.AG.boxplot.png"),
            plot=.,
            width.r=0.9, height.r=0.9
        )


        ## plot editing level (GSE133854 PG lost)
        ggplot(..to.plot.dt[completely.lost.in.GSE133854.PG.embryos == TRUE], aes(x=paste(sep="", CHROM, "_", POS), y=AF)) ->.;
        . + geom_boxplot() ->.;
        . + coord_flip() -> .;
        . + theme_pubr() ->.;
        . + theme(strip.text.y.right = element_text(angle = 0)) -> .;
        . + facet_grid(stage.description.ordered ~., scales="free_y", space="free_y") -> .;
        . + labs(x="", y="Editing level") ->.;
        . + scale_y_continuous(limits=c(0.1, 1), breaks=seq(0.1, 1, 0.3)) -> .;
        ggsave.A4(
            filename=glue("{output.directory}/all.107edits.normal.editing.level.GSE133854.PG.boxplot.png"),
            plot=.,
            width.r=0.9, height.r=0.9
        )


    }
    
}
