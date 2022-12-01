library("data.table")
library("eulerr")
library("magrittr")
library("foreach")
library("doMC")
library("iterators")
library("glue")
library("stringi")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA.intersection.all.edits/"
dir.create(output.directory, recursive=TRUE)

original.TargetScan.and.miRanda.match.dt <- fread(glue("{output.directory}/original.TargetScan.and.miRanda.match.dt.gz"))
original.TargetScan.and.miRanda.match.dt[, ID:=.I]

{
    
    original.TargetScan.and.miRanda.match.dt -> .;
    .[, list(predicted.by.TargetScan, predicted.by.miRanda)] -> .;
    .[is.na(predicted.by.TargetScan) == TRUE, predicted.by.TargetScan := FALSE] -> .;
    .[is.na(predicted.by.miRanda) == TRUE, predicted.by.miRanda := FALSE] -> .;
    .[, list("Predicted by\nTargetScan"=predicted.by.TargetScan, "Seed region\npredicted by miRanda"=predicted.by.miRanda)] -> .;
    ##
    . -> ..to.plot.dt
    ##
    png(glue("{output.directory}/original.TargetScan.and.miRanda.venn.png"), width=12, height=8, units="cm", res=600)
    print(plot(euler(..to.plot.dt), quantities=list(fontsize=7), legend=list(labels=c("Predicted by\nTargetScan", "Seed region\npredicted by miRanda")), fills=c("lightblue", "#FF8F8F")))
    dev.off()   
    
}
