library("data.table")
library("ggVennDiagram")
library("readxl")
library("foreach")
library("glue")
source("./scripts/common/ggpubr.A4.R")

input.edits.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")



## within-normal overlap, and within-abnormal overlap, boxplot
{


    ## 1.1. get valid sample groups
    {
        input.edits.dt -> .;
        .[stage %in% c('zygote', 'zygote.2PN', '2-cell', '4-cell', '8-cell', 'morula')] -> .;
        ## get sample info
        unique(.[, list(SAMPLE, stage, is.normal, treatment, disease)]) -> .;
        ## keep only those sample groups with >=3 samples
        .[, list(sample.count=.N), list(stage, is.normal, treatment, disease)][sample.count >= 3] -> .;
        .
    } -> ..valid.sample.group.dt

    ## 1.2. get details of these sample groups
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
        ## prettify the labels
        .[, stage.prettified:=stage][stage=='zygote.2PN', stage.prettified:="zygote (2PN)"] -> .;
        .[, disease.prettified:=c("androgenetic"="/AG", "parthenogenetic"="/PG", "viability predicted to be good"="/good viab.", "viability predicted to be bad"="/bad viab.")[disease]][is.na(disease.prettified) == TRUE, disease.prettified:=""] -> .;
        .[, stage.and.is.normal.and.disease.prettified:=paste(sep="", stage.prettified, "/", c("abnormal", "normal")[is.normal+1], disease.prettified)][, stage.and.is.normal.and.disease.prettified.ordered:=factor(stage.and.is.normal.and.disease.prettified, levels=c("zygote/normal", "zygote/abnormal/AG", "zygote/abnormal/PG", "2-cell/normal", "2-cell/abnormal/AG", "2-cell/abnormal/PG", "4-cell/normal", "4-cell/abnormal/AG", "4-cell/abnormal/PG", "8-cell/normal", "8-cell/abnormal/AG", "8-cell/abnormal/PG", "morula/normal", "morula/abnormal/AG", "morula/abnormal/PG", "zygote (2PN)/normal/good viab.", "zygote (2PN)/abnormal/bad viab."))] -> .;
        .
    } -> SAMPLE.info.dt


    rbindlist(lapply(1:10, FUN=function(temp.seed){
        set.seed(temp.seed)
        SAMPLE.info.dt[, .SD[, list(SAMPLE)][sample(x=.N, size=floor(0.5*.N), replace=FALSE), split.id:="Set.1"][is.na(split.id)==TRUE, split.id:="Set.2"], list(stage, stage.and.is.normal.and.disease.prettified)] -> TEMP.dt
        TEMP.dt[, seed:=temp.seed]
    }), use.names=TRUE) -> sampling.strategy.dt
           
    rbindlist(foreach(temp.seed=1:10) %do% {
        print(glue("{date()} : processing temp seed {temp.seed}"))
        sampling.strategy.dt[seed==temp.seed] -> TEMP.sampling.strategy.dt
        merge(x=input.edits.dt, y=TEMP.sampling.strategy.dt[, list(SAMPLE, stage.and.is.normal.and.disease.prettified, split.id)], by=c("SAMPLE"), all.x=FALSE, all.y=TRUE) -> TEMP.dt
        unique(TEMP.dt[, list(ID, stage, stage.and.is.normal.and.disease.prettified, split.id, hit=TRUE)]) -> TEMP.dt
        dcast(TEMP.dt, stage + stage.and.is.normal.and.disease.prettified + ID ~ split.id, value.var="hit", fill=FALSE) -> TEMP.dt
        TEMP.dt[, final.set:=c("Set.1.only", "Set.2.only", "both")[Set.1 * 1 + Set.2 * 2]] -> TEMP.dt
        TEMP.dt[, list(count=.N), list(stage, stage.and.is.normal.and.disease.prettified, final.set)] -> TEMP.dt
        TEMP.dt[, `:=`(pct=count/sum(count) * 100), list(stage, stage.and.is.normal.and.disease.prettified)] -> TEMP.dt
        TEMP.dt[, seed:=temp.seed] -> TEMP.dt
        TEMP.dt
    }, use.names=TRUE) -> all.overlap.pcts.dt

    all.overlap.pcts.dt[final.set=="both"] -> to.plot.dt
    
    ggplot(to.plot.dt, aes(x=stage.and.is.normal.and.disease.prettified, y=pct)) -> .;
    . + geom_boxplot() -> .;
    . + coord_flip() -> .;
    . + theme_pubr() -> .;
    . + lims(y=c(0, 50)) -> .;
    . + labs(x="", y="Percentage of overlap\nbetween the two within-group subsets (%)") -> .;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/subset.overlap.santiy.check.png",
        plot=.,
        width.r=0.9, height.r=0.6
    )    
}

