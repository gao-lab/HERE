library("data.table")
library("readxl")
library("glue")
library("missForest")
library("umap")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/F2A.dimensional.reduction/"
dir.create(output.directory, recursive=TRUE)

input.edits.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")
setkey(input.edits.dt, SAMPLE)


## total count per stage, venn diagram
{

    ## 1.2. get edit info on valid samples
    {
        input.edits.dt -> .;
        ## keep stages of interest
        .[stage %in% c('zygote', 'zygote.2PN', '2-cell', '4-cell', '8-cell', 'morula')] -> .;
        ## keep samples related to F2A depth examination
        .[(is.normal == TRUE) | (disease %in% c("viability predicted to be bad", "androgenetic", "parthenogenetic"))] -> .;
        ##
        .
    } -> ..edit.profile.of.valid.samples.to.reduce.dimension.dt
    fwrite(..edit.profile.of.valid.samples.to.reduce.dimension.dt, glue("{output.directory}/edit.profile.of.valid.samples.to.reduce.dimension.dt.gz"))
    ##
    ..edit.profile.of.valid.samples.to.reduce.dimension.dt
} -> edit.profile.of.valid.samples.to.reduce.dimension.dt


## umap strategy
{

    edit.profile.of.valid.samples.to.reduce.dimension.dt -> .;
    dcast(data.table(., constant=1), SAMPLE ~ ID, value.var="constant", fill=0) -> .;
    date(); umap(.[, -1], verbose=TRUE) -> ..umap.result ; date()
    
}

## GTEx editome strategy
{

    copy(edit.profile.of.valid.samples.to.reduce.dimension.dt) -> .;
    .[, `:=`(count.of.ge300.ANs=sum(AN>=300), count.of.samples=.N), list(ID)] -> .;
    .[, percentage.of.samples.with.ge300.AN:=count.of.ge300.ANs/count.of.samples] -> .;
    ## filter for deeply sequenced editing sites
    .[AN>=300] -> .;
    ## filter for editing sites discovered in many samples
    .[percentage.of.samples.with.ge300.AN>=1/3] -> .;
    ##
    . -> ..filtered.dt
    
    ## reshape into wide format and impute the editing sites
    dcast(..filtered.dt, SAMPLE ~ ID, value.var="AF", fill=0) -> ..filtered.wide.dt;
    as.matrix(..filtered.wide.dt[, -1]) -> .;
    ..filtered.wide.dt[, SAMPLE] -> rownames(.);
    missForest(.)$ximp -> ..filtered.wide.imputed.matrix
    saveRDS(..filtered.wide.imputed.matrix, glue("{output.directory}/filtered.wide.imputed.matrix.RDS"))

    ## prcomp
    prcomp(..filtered.wide.imputed.matrix) -> .
    data.table(SAMPLE=rownames(.$x), .$x[, c("PC1", "PC2")]) -> .;
    merge(
        x=.,
        y=unique(edit.profile.of.valid.samples.to.reduce.dimension.dt[, list(SAMPLE, normal.type=c("others", "normal")[is.normal + 1], stage)]),
        by="SAMPLE",
        all.x=TRUE, all.y=FALSE) -> .;
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ##
    . -> ..to.plot.prcomp.dt
    ggplot(..to.plot.prcomp.dt, aes(x=PC1, y=PC2, color=stage.description.ordered)) -> .;
    . + geom_point() -> .;
    . + facet_grid(~normal.type) -> .;
    . + theme_pubr() -> .;
    . + ggtitle("PCA plot") -> .;
    . + labs(color="") -> .;
    ggsave.A4(filename=glue("{output.directory}/sample.DR.scatterplot.prcomp.png"), plot=., width.r=0.75, height.r=0.5)

    ## umap
    umap(..filtered.wide.imputed.matrix, verbose=TRUE) -> .
    data.table(SAMPLE=rownames(.$layout), .$layout) -> .;
    setnames(., c("SAMPLE", "PC1", "PC2")) -> .;
    merge(
        x=.,
        y=unique(edit.profile.of.valid.samples.to.reduce.dimension.dt[, list(SAMPLE, normal.type=c("others", "normal")[is.normal + 1], stage)]),
        by="SAMPLE",
        all.x=TRUE, all.y=FALSE) -> .;
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ##
    . -> ..to.plot.umap.dt
    ggplot(..to.plot.umap.dt, aes(x=PC1, y=PC2, color=stage.description.ordered)) -> .;
    . + geom_point() -> .;
    . + facet_grid(~normal.type) -> .;
    . + theme_pubr() -> .;
    . + ggtitle("UMAP plot") -> .;
    . + labs(color="") -> .;
    . + scale_color_manual(values=c("#d02c07", "#f96c4f", "#f1b785", "#82b740", "#36bdc0", "#3e78b9")) -> .;
    ggsave.A4(filename=glue("{output.directory}/sample.DR.scatterplot.umap.png"), plot=., width.r=0.75, height.r=0.5)
    


}
