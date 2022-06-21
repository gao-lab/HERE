library("data.table")
library("readxl")
library("foreach")
library("iterators")
library("magrittr")
library("clusterProfiler")
library("org.Hs.eg.db")
source("./scripts/common/ggpubr.A4.R")


RE.targeted.genes.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.gene/RE.targeted.genes.dt.csv.gz")

temp.GO.directory <- "report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.gene/GO"
dir.create(paste(sep="", temp.GO.directory), recursive=TRUE)

GO.raw.result.dt <- {
    RE.targeted.genes.dt ->.;
    ##s pick stages
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', "2-cell", '4-cell')] ->.;
    ##s pick Annotations
    .[Annotation.class.pasted.description.mixed %in% c("exonic", "intronic")] ->.;
    ## prepare for running GO
    setkey(., stage, Annotation.class.pasted.description.mixed)
    .[, Gene_ID.no.version:=sub(pattern="\\.[0-9]*$", replacement="", x=Gene_ID)]
    temp.combinations.dt <- .[, list(stage, Annotation.class.pasted.description.mixed)] %>% unique %>% as.data.frame(stringsAsFactors=FALSE) %>% merge(y=data.frame(temp.ontology=c("BP", "MF", "CC"), stringsAsFactors=FALSE)) %>% data.table
    ## run GO for each stage x Annotation
    foreach(temp.keys=iter(temp.combinations.dt, by='row')) %do%
    {
        cat(date(), " Processing ", temp.keys[1, stage], " x ", temp.keys[1, Annotation.class.pasted.description.mixed], " x ", temp.keys[1, temp.ontology], "\n")
        temp.path <- paste(sep="", temp.GO.directory, "/", temp.keys[1, stage], "_", temp.keys[1, Annotation.class.pasted.description.mixed], "_", temp.keys[1, temp.ontology], "_enrichGO.result.RDS")
        if (file.info(temp.path)$size %in% c(NA, 0)==TRUE){
            temp.enrichGO.result <- enrichGO(.[temp.keys[, list(stage, Annotation.class.pasted.description.mixed)], Gene_ID.no.version], OrgDb='org.Hs.eg.db', keyType="ENSEMBL", ont=temp.keys[1, temp.ontology], pvalueCutoff=0.01)
            print(temp.enrichGO.result)
            saveRDS(temp.enrichGO.result, temp.path)
        }
        temp.enrichGO.result <- readRDS(temp.path)
        if (is.null(temp.enrichGO.result) == TRUE){
            data.table(temp.keys)
        } else {
            temp.enrichGO.result@result %>% data.table(temp.keys) %>% {.[pvalue<0.1]} ## increases statistical power
        }
    } -> .; rbindlist(., use.names=TRUE, fill=TRUE) ->.;
}


fwrite(GO.raw.result.dt, "report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.gene/GO.raw.result.dt.csv.gz")


{
    copy(GO.raw.result.dt) -> .;
    ##= adjust p-values
    .[, p.adjusted.BH:=p.adjust(p=pvalue, method="BH")]
    ##s keep those combinations without enriched terms for double checking
    .[p.adjusted.BH<0.1 | is.na(p.adjusted.BH)] -> .;
    ##s pick statistically significant results
    .[p.adjusted.BH<0.05] ->.;
    ##c count occurrence of GO terms
    .[, Description.occurrence:=.N, list(Description, Annotation.class.pasted.description.mixed)][Description.occurrence>=3] -> .;
    ## prettify Description
    .[, Description.renamed:=Description]
    .[Description == "DNA-dependent DNA replication maintenance of fidelity", Description.renamed:="DNA-dependent DNA\nreplication maintenance of fidelity"]
    .[Description == "double-strand break repair via homologous recombination", Description.renamed:="double-strand break repair\nvia homologous recombination"]
    .[Description == "regulation of viral genome replication", Description.renamed:="regulation of\nviral genome replication"]
    .[Description == "negative regulation of DNA replication", Description.renamed:="negative regulation of\nDNA replication"]
    ## prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, before.or.after.8.cell.ordered:=factor(before.or.after.8.cell, levels=unique(temp.stage.dt[, before.or.after.8.cell]))] ->.;
    ## prettify p-values
    .[, p.adjusted.BH.cut:=cut(x=p.adjusted.BH, breaks=c(0, 0.01, 0.05, 0.1, 1))]
    ## prettify Annotation
    .[, Annotation.class.pasted.description.mixed.renamed:=c("exonic"="always exonic\n(mostly 3'-UTR)", "intronic"="always intronic")[Annotation.class.pasted.description.mixed]]
    ## start plotting
    ggplot(., aes(x=stage.description.ordered, y = Description.renamed, fill=p.adjusted.BH.cut)) ->.;
    . + geom_tile() ->.;
    . + facet_grid(temp.ontology ~ Annotation.class.pasted.description.mixed.renamed, scale="free_y", space="free") ->.;
    . + labs(x="", y="", fill="BH-adjusted\np-value") -> .;
    . + theme_pubr(base_size=10) ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(lineheight=1)) ->.;
    . + guides(fill=guide_legend(nrow=2)) ->.;
    . + scale_fill_manual(values=c("#F8766D", "#AB514B")) ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.gene/RE.gene.GO.tileplot.png",
        plot=.,
        width.r=0.45, height.r=0.6
    )    
}
