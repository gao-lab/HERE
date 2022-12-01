library("data.table")
library("ggVennDiagram")
library("readxl")
source("./scripts/common/ggpubr.A4.R")

input.edits.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")

input.variants.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz")

## total count, venn diagram
{
    input.edits.dt -> .;
    ## generate input for ggVennDiagram
    list(
        "Normal samples"=.[is.normal==TRUE, ID],
        "Abnormal samples"=.[is.normal==FALSE, ID]
    ) ->.;
    ##
    ggVennDiagram(x=., edge_size=0.5) ->.;
    . + scale_fill_gradient(low="white", high="white") ->.;
    . + scale_color_manual(values=c("black", "black")) ->.;
    . + theme(legend.position="none") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/all.edits.venn.resized.png",
        plot=.,
        width.r=0.45, height.r=0.18
    )    
}


## per-stage count, barplot
{
    input.edits.dt ->.
    ##c get per-group edit counts
    unique(.[, list(ID, stage)]) ->.;
    .[, list(count=.N), list(stage)] -> .;
    ## prettify stage
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=rev(temp.stage.dt[, stage.description]))] ->.;
    ##plot count of edits per group
    ggplot(., aes(x=stage.description.ordered, y=count)) ->.;
    . + geom_bar(stat="identity") ->.;
    . + scale_y_log10() ->.;
    . + coord_flip() ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1)) ->.;
    . + labs(x="", y="# edits\nidentified") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/all.edits.count.per.stage.barplot.png",
        plot=.,
        width.r=0.45, height.r=0.7
    )       
}

## A-to-G ratio, boxplot
{
    
    input.variants.dt ->.
    ##c get variant count per sample x variant type
    .[, list(count=.N), list(SAMPLE, stage, event.summary)] ->.;
    ##= get variant ratio per sample
    .[, total.count:=sum(count), list(SAMPLE, stage)]
    .[, percentage:=count/total.count]
    .[, percentage.to.plot:=percentage*100]
    ##s get A-to-G variants only
    .[event.summary %in% c('A>G', 'A>G;T>C')] ->.;
    ##c get A-to-G ratio
    .[, list(A.to.G.ratio=sum(percentage.to.plot)), list(SAMPLE, stage)] ->.;
    ## add category
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))] ->.;
    . -> to.plot.dt;
    fwrite(to.plot.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/source.data.for.all.edits.A.to.G.ratio.boxplot.resized.png.csv.gz")

    ##plot A-to-G ratio
    ggplot(., aes(x="a", y=A.to.G.ratio, color=category.ordered)) ->.;
    . + geom_boxplot() ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right") ->.;
    . + guides(color=guide_legend(nrow=3)) ->.;
    . + scale_y_continuous(limits=c(0, 100)) ->.;
    . + labs(x="", y="% A-to-G\nvariants", color="") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/all.edits.A.to.G.ratio.boxplot.resized.png",
        plot=.,
        width.r=0.45, height.r=0.12
    )       
}



## Alu ratio, boxplot
{
    
    input.edits.dt ->.
    ##c get edit count per sample x SUBSET
    .[, list(count=.N), list(SAMPLE, stage, SUBSET)] ->.;
    ##= get variant ratio per sample
    .[, total.count:=sum(count), list(SAMPLE)]
    .[, percentage:=count/total.count]
    .[, percentage.to.plot:=percentage*100]
    ##s get Alu records only
    .[SUBSET == "Alu"] ->.;
    ## add category
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))] ->.;
    . -> to.plot.dt;
    fwrite(to.plot.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/source.data.for.all.edits.Alu.ratio.boxplot.resized.png.csv.gz")

    to.plot.dt -> .;
    ##plot Alu ratio
    ggplot(., aes(x="a", y=percentage.to.plot, color=category.ordered)) ->.;
    . + geom_boxplot() ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right") ->.;
    . + guides(color=guide_legend(nrow=3)) ->.;
    . + scale_y_continuous(limits=c(0, 100)) ->.;
    . + labs(x="", y="% edits\non Alu", color="") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/all.edits.Alu.ratio.boxplot.resized.png",
        plot=.,
        width.r=0.45, height.r=0.12
    )       
}
