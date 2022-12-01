library("data.table")
library("readxl")
library("foreach")
library("glue")
library("scales")
source("./scripts/common/ggpubr.A4.R")


output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.and.expression.using.miRNA.intersection.all.edits"
median.FPKM.per.gene.and.stage.dt <- fread(glue("{output.directory}/median.FPKM.per.gene.and.stage.dt.gz"))
median.AF.per.gene.and.stage.dt <- fread(glue("{output.directory}/median.AF.per.gene.and.stage.dt.gz"))

stage.and.subsequent.stages.list <- list(
    'oocyte.GV'=c('oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell'),
    'oocyte.MII'=c('zygote', '2-cell', '4-cell', '8-cell'),
    'zygote'=c('2-cell', '4-cell', '8-cell'),
    '2-cell'=c('4-cell', '8-cell'),
    '4-cell'=c('8-cell')
)


{
    
    copy(median.FPKM.per.gene.and.stage.dt) -> .;
    ## create unique label
    unique(median.AF.per.gene.and.stage.dt[, list(stage, Gene_ID, has.REE.in.this.stage=TRUE)]) -> ..has.REE.dt
    ## add REE existence
    merge(x=., y=..has.REE.dt, by=c("stage", "Gene_ID"), all.x=TRUE, all.y=FALSE) -> .;
    .[is.na(has.REE.in.this.stage) == TRUE, has.REE.in.this.stage := FALSE] -> .;
    ## prettify stage
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, before.or.after.8.cell.ordered:=factor(before.or.after.8.cell, levels=unique(temp.stage.dt[, before.or.after.8.cell]))] ->.;
    ##
    . -> ..to.plot.dt

    ## plot overall boxplot trends
    foreach(temp.stage=c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell')) %do% {
        ##
        print(glue("{date()} : plotting for {temp.stage} ..."))
        copy(..to.plot.dt) -> .;
        ## add oocyte.GV has.REEs
        .[stage==temp.stage, list(Gene_ID, has.REE.in.the.starting.stage=factor(c("Does not have REE\nin the starting stage", "Has REE in the\nstarting stage")[has.REE.in.this.stage + 1], levels=c("Has REE in the\nstarting stage", "Does not have REE\nin the starting stage")))] -> ..temp.dt
        merge(x=., y=..temp.dt, by="Gene_ID", all.x=TRUE, all.y=FALSE) -> .;
        .[Gene_ID %in% ..to.plot.dt[stage==temp.stage & FPKM.median > 0.1, Gene_ID]] -> .
        ## label the starting and succeeding stages
        .[stage==temp.stage, stage.label:="starting stage"] -> .;
        .[stage %in% stage.and.subsequent.stages.list[[temp.stage]] , stage.label:="succeeding stage"]
        .[is.na(stage.label) == TRUE, stage.label:="other"] -> .;
        ## 
        . -> ..to.subplot.dt
        ..to.subplot.dt[, list(boxplot.median=median(FPKM.median)), list(stage.description.ordered, has.REE.in.the.starting.stage, stage.label)] -> ..boxplot.median.dt
        ##
        ggplot(..to.subplot.dt, aes(x=stage.description.ordered, y=log10(FPKM.median + 1e-4), fill=stage.label)) -> .;
        . + geom_boxplot() -> .;
        . + geom_text(data=..boxplot.median.dt, mapping=aes(x=stage.description.ordered, y=-5, label=round(boxplot.median, 2))) -> .;
        . + facet_grid(has.REE.in.the.starting.stage~.) -> .;
        . + theme_pubr(base_size=10) ->.;
        . + labs(x="", fill="Temporal stage annotation\nw.r.t. the observation of REE") -> .;
        . + scale_y_continuous(limits=c(-5.2, 4.5)) -> .;
        if (temp.stage == 'oocyte.GV') {
            . + scale_fill_manual(values=c("Green", "SkyBlue"))-> .;
        } else {
            . + scale_fill_manual(values=c("Gray", "Green", "SkyBlue"))-> .;
        }
        ggsave.A4(
            filename=glue("{output.directory}/developmental.median.FPKM.on.{temp.stage}.REE.existence.png"),
            plot=.,
            width.r=0.8, height.r=0.3
        )
        ##       
    } -> not.used.variable

    ## plot heatmap of oocyte.GV-highly expressed genes
    copy(..to.plot.dt) -> .;
    'oocyte.GV' -> temp.stage
    'oocyte (GV)' -> temp.stage.description
    .[Gene_ID %in% .[stage==temp.stage & FPKM.median > 200, Gene_ID]] -> .;
    .[stage==temp.stage, list(Gene_ID, has.REE.in.the.starting.stage=factor(c(glue("Does not have\nREE in the stage\n{temp.stage.description}"), glue("Has REE in\nthe stage\n{temp.stage.description}"))[has.REE.in.this.stage + 1], levels=c(glue("Has REE in\nthe stage\n{temp.stage.description}"), glue("Does not have\nREE in the stage\n{temp.stage.description}"))))] -> ..temp.dt
    merge(x=., y=..temp.dt, by="Gene_ID", all.x=TRUE, all.y=FALSE) -> .;
    . -> to.plot.heatmap.dt;
    ##
    ## order genes by hclust
    {
        dcast(to.plot.heatmap.dt, Gene_ID ~ stage, value.var="FPKM.median") -> temp.dt;
        data.frame(temp.dt[, -1]) -> .;
        temp.dt[, Gene_ID] -> rownames(.);
        hclust(dist(.)) -> .;
        temp.dt[.$order, Gene_ID] -> Gene_ID.levels
    }
    to.plot.heatmap.dt[, Gene_ID.ordered:=factor(Gene_ID, levels=Gene_ID.levels)]
    ##
    ## add lm coef and pval.adj
    {
        to.plot.heatmap.dt -> .;
        dcast(., has.REE.in.the.starting.stage + Gene_ID.ordered ~ stage, value.var="FPKM.median") -> .;
        .[, list(has.REE.in.the.starting.stage, Gene_ID.ordered, oocyte.GV, oocyte.MII, zygote, `2-cell`, `4-cell`, `8-cell`)] -> TEMP.raw.dt;
        {
            TEMP.raw.dt[,  list(oocyte.GV, oocyte.MII, zygote, `2-cell`, `4-cell`, `8-cell`)] -> .;
            data.table(do.call(rbind, apply(., MARGIN=1, FUN=function(temp.y){1:6 -> temp.x; lm(temp.y~temp.x) -> temp.lm; coef(summary(temp.lm)) -> temp.coef; temp.coef["temp.x", c("Estimate", "Pr(>|t|)")]}, simplify=FALSE))) -> .;
            . -> TEMP.lm.dt
        }
        data.table(TEMP.raw.dt[, list(has.REE.in.the.starting.stage, Gene_ID.ordered)], TEMP.lm.dt) -> .;
        .[, `sign of\nregression\ncoefficient`:=c("-1"="negative", "1"="positive")[as.character(sign(Estimate))]] -> .;
        .[, pval.adj:=p.adjust(`Pr(>|t|)`, method="BH")] -> .;
        .[, `adjusted\np-value`:=c(">=0.1", "<0.1")[(pval.adj<0.1) + 1]]
        melt(., id.vars=c("has.REE.in.the.starting.stage", "Gene_ID.ordered"), measure.vars=c("sign of\nregression\ncoefficient", "adjusted\np-value"), variable.name="metric", value.name="value") -> .;
        .[, metric.ordered:=factor(metric, levels=c("adjusted\np-value", "sign of\nregression\ncoefficient"))] -> .;
        .[, value.ordered:=factor(value, levels=c(">=0.1", "<0.1", "positive", "negative"))]
        . -> to.plot.lm.heatmap.dt
    }
    ##
    ## plot expression heatmap
    ggplot(to.plot.heatmap.dt, aes(x=stage.description.ordered, y=Gene_ID.ordered, fill=log10(FPKM.median + 1e-4))) -> .;
    . + geom_tile() -> .;
    . + facet_grid(has.REE.in.the.starting.stage~., scales="free_y", space="free_y") -> .;
    . + theme_pubr(base_size=10) ->.;
    . + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text.y=element_blank()) -> .;
    . + labs(x="", y="Gene with FPKM median > 200 in oocyte (GV)") ->.;
    . + scale_fill_gradientn(colours=c("blue", "lightblue", "white", "#FF8888", "red"), values=rescale(c(-4,  0, 0.2, 3,5))) -> .;
    . -> expression.heatmap.ggplot
    ##
    ## plot lm heatmap
    ggplot(to.plot.lm.heatmap.dt, aes(x=metric.ordered, y=Gene_ID.ordered, fill=value.ordered)) -> .;
    . + geom_tile() -> .;
    . + facet_grid(has.REE.in.the.starting.stage~., scales="free_y", space="free_y") -> .;
    . + theme_pubr(base_size=10) ->.;
    . + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y = element_blank()) -> .;
    . + labs(x="", y="", fill="") ->.;
    . + guides(fill=guide_legend(nrow=2)) -> .;
    . + scale_fill_manual(values=c("gray", "lightgreen", "#FF8888", "lightblue")) -> .;
    . -> lm.heatmap.ggplot
    ##
    ggarrange(plotlist=list(expression.heatmap.ggplot, lm.heatmap.ggplot), nrow=1, align="h", widths=c(2.5,1)) -> .;
    ggsave.A4(
        filename=glue("{output.directory}/developmental.median.FPKM.on.oocyte.GV.heatmap.png"),
        plot=.,
        width.r=0.9, height.r=0.8
    )
    
}
