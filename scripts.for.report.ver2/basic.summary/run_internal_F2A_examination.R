library("data.table")
library("ggVennDiagram")
library("readxl")
library("glue")
library("foreach")
library("iterators")
library("magrittr")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/F2A.examination/"
dir.create(output.directory, recursive=TRUE)

input.edits.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")
setkey(input.edits.dt, SAMPLE)


input.variants.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz")

## total count per stage, venn diagram
{

    ## 1.1. get edit count ~ sequencing depth (srr.total.bases)
    {
        input.edits.dt -> .;
        .[stage %in% c('zygote', 'zygote.2PN', '2-cell', '4-cell', '8-cell', 'morula')] -> .;
        .[, list(edit.count=length(ID)), list(SAMPLE, gse, stage, is.normal, srr.total.bases, treatment, disease)] -> .;
    } -> ..edit.count.per.sample.info.dt
    fwrite(..edit.count.per.sample.info.dt, glue("{output.directory}/edit.count.per.sample.info.dt.gz"))

    ## 1.2. get valid sample groups
    {
        input.edits.dt -> .;
        .[stage %in% c('zygote', 'zygote.2PN', '2-cell', '4-cell', '8-cell', 'morula')] -> .;
        ## get sample info
        unique(.[, list(SAMPLE, stage, is.normal, treatment, disease)]) -> .;
        ## keep only those sample groups with >=3 samples
        .[, list(sample.count=.N), list(stage, is.normal, treatment, disease)][sample.count >= 3] -> .;
        .
    } -> ..valid.sample.group.dt

    ## 1.3. get valid sample pairs per group
    {
        ##        
        input.edits.dt -> .;
        .[stage %in% c('zygote', 'zygote.2PN', '2-cell', '4-cell', '8-cell', 'morula')] -> .;
        ## get sample info
        unique(.[, list(SAMPLE, stage, is.normal, treatment, disease)]) -> .;
        ##
        ## keep valid sample groups only
        merge(x=., y=..valid.sample.group.dt,
              by=c("stage", "is.normal", "treatment", "disease"),
              all.x=FALSE, all.y=TRUE) -> .;
        ##
        ## get sample pairs per group
        .[, setnames(data.table(t(combn(SAMPLE, 2))), c("SAMPLE.1", "SAMPLE.2")), list(stage, is.normal, treatment, disease)] -> .;
        ##
        ## append to each sample its sequencing depth
        merge(x=., y=..edit.count.per.sample.info.dt[, list(SAMPLE, srr.total.bases)],
              by.x="SAMPLE.1", by.y="SAMPLE",
              all.x=TRUE, all.y=FALSE) -> .;
        merge(x=., y=..edit.count.per.sample.info.dt[, list(SAMPLE, srr.total.bases)],
              by.x="SAMPLE.2", by.y="SAMPLE",
              all.x=TRUE, all.y=FALSE) -> .;
        ##
        .
    } -> ..valid.sample.pairs.per.group.dt

    ## 1.4. get edit overlap between each pair of samples from the same group
    {
        ##
        ..valid.sample.pairs.per.group.dt -> .
        foreach(ROW=iter(., by="row")) %do% {
            input.edits.dt[ROW[1, SAMPLE.1], ID] -> ..temp.edits.in.sample.1.vector
            input.edits.dt[ROW[1, SAMPLE.2], ID] -> ..temp.edits.in.sample.2.vector
            data.table(
                ROW,
                intersection.size=length(intersect(..temp.edits.in.sample.1.vector, ..temp.edits.in.sample.2.vector)),
                union.size=length(union(..temp.edits.in.sample.1.vector, ..temp.edits.in.sample.2.vector))
            ) -> RESULT
            RESULT
        } %>% rbindlist(use.names=TRUE) -> .
        ## compute overlap ratio
        .[, overlap.ratio:=intersection.size / union.size] -> .;
        ##
        .
    } -> ..edit.overlap.between.sample.pairs.per.group.dt
    fwrite(..edit.overlap.between.sample.pairs.per.group.dt, glue("{output.directory}/edit.overlap.between.sample.pairs.per.group.dt.gz"))
    ##
}

edit.count.per.sample.info.dt <- fread(glue("{output.directory}/edit.count.per.sample.info.dt.gz"))
edit.overlap.between.sample.pairs.per.group.dt <- fread(glue("{output.directory}/edit.overlap.between.sample.pairs.per.group.dt.gz"))

{

    ## 2.1. edit count x sequencing depth
    {
        copy(edit.count.per.sample.info.dt) -> .;
        .[, srr.total.bases.divided.by.1e9 := srr.total.bases / 1e9] -> .;
        ## prettify stages
        temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
        merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
        ## prettify normal.type
        .[, normal.type:=c("abnormal", "normal")[as.integer(is.normal) + 1]] -> .;
        ##
        . -> ..to.plot.dt;
        ggplot(..to.plot.dt, aes(x=srr.total.bases.divided.by.1e9, y=edit.count)) -> .;
        . + geom_point() -> .;
        . + facet_wrap(~stage.description.ordered + normal.type, scales="free") -> .;
        . + theme_pubr() -> .;
        . + theme(axis.text.x=element_text(angle=45, hjust=1)) -> .;
        . + labs(x="# total bases sequenced / 1E9", y="# edits identified") ->.;       
        ggsave.A4(filename=glue("{output.directory}/edit.counts.vs.sequencing.depth.scatterplot.png"), plot=., width.r=0.75, height.r=0.5)
    }

    ## 2.2. edit count x sequencing depth
    {
        
        copy(edit.overlap.between.sample.pairs.per.group.dt) -> .;
        .[, log10.srr.total.bases.distance := abs(log10(srr.total.bases.x) - log10(srr.total.bases.y))] -> .;
        .[, stage.prettified:=stage][stage=='zygote.2PN', stage.prettified:="zygote (2PN)"] -> .;
        .[, spearman.cor:=cor(log10.srr.total.bases.distance, overlap.ratio, method="spearman"), list(stage, is.normal, disease)][, spearman.cor.rounded:=paste(sep="", "Spearman cor.: ", round(spearman.cor, 4))] -> .;
        .[, disease.prettified:=c("androgenetic"="/AG", "parthenogenetic"="/PG", "viability predicted to be good"="/good viab.", "viability predicted to be bad"="/bad viab.")[disease]][is.na(disease.prettified) == TRUE, disease.prettified:=""] -> .;
        .[, stage.and.is.normal.and.disease.prettified:=paste(sep="", stage.prettified, "/", c("abnormal", "normal")[is.normal+1], disease.prettified)][, stage.and.is.normal.and.disease.prettified.ordered:=factor(stage.and.is.normal.and.disease.prettified, levels=c("zygote/normal", "zygote/abnormal/AG", "zygote/abnormal/PG", "2-cell/normal", "2-cell/abnormal/AG", "2-cell/abnormal/PG", "4-cell/normal", "4-cell/abnormal/AG", "4-cell/abnormal/PG", "8-cell/normal", "8-cell/abnormal/AG", "8-cell/abnormal/PG", "morula/normal", "morula/abnormal/AG", "morula/abnormal/PG", "zygote (2PN)/normal/good viab.", "zygote (2PN)/abnormal/bad viab."))] -> .;
        . -> ..to.plot.dt;
        ##
        ggplot(..to.plot.dt, aes(x=log10.srr.total.bases.distance, y=overlap.ratio)) -> .;
        . + geom_point() -> .;
        . + geom_smooth(formula = y~x, color="orange") -> .;
        . + facet_wrap(~stage.and.is.normal.and.disease.prettified.ordered + spearman.cor.rounded, scales="free", ncol=3) -> .;
        . + theme_pubr() -> .;
        . + theme(axis.text.x=element_text(angle=45, hjust=1)) -> .;
        . + labs(x="|log10(# total bases sequenced of A) - log10(# total bases sequenced of B)|\nfor a given sample pair (A, B)", y="overlap ratio (= (# edits identified in both samples) / (# edits identified in A, or B, or both) )\nfor a given sample pair (A, B)") ->.;       
        ggsave.A4(filename=glue("{output.directory}/edit.overlap.vs.sequencing.depth.distance.scatterplot.png"), plot=., width.r=0.98, height.r=0.95)
        
    }

}
