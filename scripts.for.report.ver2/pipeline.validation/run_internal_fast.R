library("data.table")
library("ggpubr")
library("foreach")
library("iterators")
source("./scripts/common/ggpubr.A4.R")

fread("./report.ver2/pipeline.validation/210203-GSE144296.A375/combined.RNA.DNA.comparison.dt.csv.gz") -> combined.RNA.DNA.comparison.dt


{
    
    combined.RNA.DNA.comparison.dt ->.;
    ## get counts of variants per sample
    .[, list(count=.N), list(SAMPLE, is.DNA.detected, RNA.edit.type)] ->.;
    ## add all missing combinations
    setkey(., SAMPLE, is.DNA.detected, RNA.edit.type)
    .[CJ(SAMPLE, is.DNA.detected, RNA.edit.type, unique=TRUE)] ->.;
    setnafill(x=., type="const", fill=0, cols="count")
    ## prettify is.DNA.detected
    .[, DNA.support:=c("Not overlapping with genomic variants", "Overlapping with genomic variants")[(is.DNA.detected == TRUE) + 1]] -> .
    ## take the mean across samples
    .[, list(mean.count=mean(count)), list(DNA.support, RNA.edit.type)] -> .;
    ##
    . -> ..to.plot.dt
    
    ## start plotting
    ggplot(..to.plot.dt, aes(x="a", y=mean.count, fill=DNA.support)) ->.;
    ## plot boxplots
    . + geom_bar(stat="identity", position="dodge") ->.;
    ## facet
    . + facet_grid(RNA.edit.type~., scales="free") ->.;
    ## add theme
    . + theme_pubr(base_size=10) ->.;
    . + theme(axis.text.x=element_blank()) ->.;
    . + guides(fill=guide_legend(nrow=2)) -> .;
    . + labs(x="A>G and A>G;T>C variants", y="Mean(# variants) across samples", fill="") ->.;
    ## save image
    ggsave.A4(
        filename="./report.ver2/pipeline.validation/210203-GSE144296.A375/validation.resized.bar.plot.png",
        plot=.,
        width.r=0.35, height.r=0.25)
    
}
