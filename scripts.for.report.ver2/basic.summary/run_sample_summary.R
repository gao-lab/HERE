library("data.table")
library("readxl")
library("glue")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/sample.summary/"
dir.create(output.directory, recursive=TRUE)


{
    fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt") -> .;
    ## pick valid samples
    .[( gse=="GSE36552" ) | (srr.mean.avgspotlen >= 150)] -> .;
    ## compute counts
    .[, list(count=.N), list(stage, is.normal)] -> .;
    ## prettify normal type
    .[, normal.type:=c("abnormal", "normal")[as.integer(is.normal) + 1]] -> .;
    ## prettify stage
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=rev(temp.stage.dt[, stage.description]))] ->.;
    ##
    . -> ..to.plot.dt
    fwrite(..to.plot.dt[, list(stage=gsub(pattern="(\r|)\n", replacement=" ", x=stage.description.ordered), normal.type, count.of.samples=count)], glue("{output.directory}/sample.count.csv"))
    ##
    ggplot(..to.plot.dt, aes(x=stage.description.ordered, y=log10(count+0.1))) ->.;
    . + geom_bar(stat="identity") ->.;
    . + coord_flip() ->.;
    . + facet_grid(~normal.type, scales="free_x") -> .;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) ->.;
    . + labs(x="", y="log10(# samples + 0.1)") ->.;
    ggsave.A4(
        filename=glue("{output.directory}/sample.summary.barplot.png"),
        plot=.,
        width.r=0.9, height.r=0.7
    )       
    ##
}
