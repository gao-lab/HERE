#' ---
#' title: "F1D summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("scales")
library("ggtext")
source("./scripts/common/ggpubr.A4.R")

merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")



for.plot.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt %>%
    {.[POS==37519170 & CHROM=='chr20']} %>%
    {.[is.normal==TRUE & is.na(AC)==FALSE]} %>%
    {.[, list(CHROM, POS, SAMPLE, gse, stage, AF, AC, AN)]} %>%
    {.[stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell", "8-cell", "morula", "hESC", "trophoblast", "ICM")]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=c("gene", "genomic DNA", "unedited RNA", "edited RNA", "unedited AA", "edited AA", temp.stage.dt[, stage.description]))][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]}

fwrite(for.plot.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/sourec.data.for.BLCAP.example.resized.png.csv.gz")

BLCAP.example.ggplot <- for.plot.dt %>%
    {
        ggplot(., aes(x=paste(sep="", CHROM, ":", POS), y=SAMPLE, fill=cut(x=AF, breaks=c(0, 0.1, 0.2, 0.4, 1)))) +
            geom_tile() +
            labs(x="", y="", fill="editing level") +
            geom_blank(data=data.table(POS=seq(37519140, 37519200, 4), SAMPLE="GSM3928470", stage.description.ordered=factor("morula", levels=levels(.[, stage.description.ordered])), CHROM="chr20", AF=as.numeric(NA))) +
            facet_grid(stage.description.ordered~., scales="free", switch="y") +
            theme_pubr() +
            theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.background = element_rect(color="black", fill="white"), strip.text.y.left = element_text(angle = 0), plot.margin=margin(t=75,r=0,b=0,l=0), legend.position=c(0.3, 1.1)) +
            guides(fill=guide_legend(nrow=1, title.hjust=0.5)) +
            scale_x_discrete(breaks=c("chr20:37519170")) +
            geom_richtext(data=data.table(CHROM="chr20", POS=37519170, SAMPLE=0, text=c("<i>BLCAP</i>", "<span style='font-family:monospace'>--GCAATACAT-></span>", "<span style='font-family:monospace'><-CGUUAUGUA--</span>", "<span style='font-family:monospace'><-CGUU<span style='color:#ff0000'>G</span>UGUA--</span>", "<span style='font-family:monospace'><--C--Y--M---</span>", "<span style='font-family:monospace'><--C--<span style='color:#ff0000'>C</span>--M---</span>"), stage.description.ordered=factor(c("gene", "genomic DNA", "unedited RNA", "edited RNA", "unedited AA", "edited AA"), levels=levels(.[, stage.description.ordered])), AF=as.numeric(NA)), mapping=aes(x=paste(sep="", CHROM, ":", POS), y=SAMPLE, label=text), fill=NA, label.size=0) + scale_color_manual(values=c("black", "red")) + guides(color="none") 
    }
ggsave.A4(filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/BLCAP.example.resized.png", plot=BLCAP.example.ggplot, width.r=0.4, height.r=0.6)
