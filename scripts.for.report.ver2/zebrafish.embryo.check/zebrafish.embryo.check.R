library("readxl")
library("data.table")
source("./scripts/common/ggpubr.A4.R")

dir.create("report.ver2/zebrafish.embryo.check", recursive=TRUE)

{

    "./external/zebrafish.embryo/SuppTable2.xlsx" -> .;
    read_excel(.) -> .;
    as.data.table(.) -> .;
    . -> to.plot.dt
    
    ggplot(to.plot.dt, aes(x=editingMean)) -> .;
    . + geom_histogram() -> .;
    . + theme_pubr() -> .;
    . + labs(x="Mean editing level", y="Count") ->.;    
    ggsave.A4(
        filename="./report.ver2/zebrafish.embryo.check/149.coding.edits.editing.level.histogram.png",
        plot=.,
        width.r=0.45, height.r=0.2
    )    

    
    
    
}
