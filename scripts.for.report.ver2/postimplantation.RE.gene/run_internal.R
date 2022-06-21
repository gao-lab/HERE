library("data.table")
library("readxl")
library("ggalluvial")
library("foreach")
library("magrittr")
source("./scripts/common/ggpubr.A4.R")

RE.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

total.sample.count.for.normal.stages.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.sample.count.for.normal.stages.dt.csv")
total.sample.count.for.normal.stages.mapping.vector <- {
    total.sample.count.for.normal.stages.dt -> .;
    temp.vector <- .[, total.sample.count];
    names(temp.vector) <- .[, stage]
    temp.vector
}



## generate table of RE-targeted genes and plot counts per stage
postimplantation.RE.genes.and.their.RE.dt <- {
    RE.valid.genes.only.dt -> .;
    ##s pick edited records only
    .[is.na(AC)==FALSE] ->.;
    ##c collapse Annotation.class 
    .[, list(Annotation.class.pasted=paste(collapse=";", Annotation.class %>% sort %>% unique)), list(Gene_ID, Gene_Name, SAMPLE, gse, stage)] ->.;
    .[, Annotation.class.pasted.description:=c('exonic.or.splicing.related'='exonic', 'purely.intronic'='intronic', 'exonic.or.splicing.related;purely.intronic'='mixed')[Annotation.class.pasted]] ->.;
    ##= compute % occurrence of each gene x Annotation.class.pasted.
    .[, total.sample.count:=total.sample.count.for.normal.stages.mapping.vector[stage]]
    .[, gene.occurrence:=.N, list(Gene_ID, Gene_Name, stage)]
    .[, gene.occurrence.pct:=gene.occurrence/total.sample.count]
    ##c pick those with gene occurrence pct >= 80% samples. Note that the resulting `dt` conforms to alluvial form (i.e., <=1 record for each gene x stage)
    .[gene.occurrence.pct>=0.8] ->.;
    ##c collapse to each stage
    .[, list(Annotation.class.pasted.description.mixed=paste(collapse=";", Annotation.class.pasted.description %>% sort %>% unique)), list(stage, Gene_ID, Gene_Name)] ->.;
    ##= reduce Annotation to exon/intron/mixed
    .[, Annotation.class.pasted.description.mixed.ordered:=factor(c("exonic"="always exonic\n(mostly 3'-UTR)", "intronic"="always intronic", "mixed"="mixed", "exonic;intronic"="mixed", "exonic;mixed"="mixed", "intronic;mixed"="mixed", "exonic;intronic;mixed"="mixed")[Annotation.class.pasted.description.mixed], levels=c("always exonic\n(mostly 3'-UTR)", "always intronic", "mixed"))] 
    ##s only listing those 'normal' stages with sample size >= 10
    .[stage %in% c('morula', 'blastocyst.late', 'trophoblast', 'ICM', 'CTB', 'STB', 'MTB', 'EVT', 'epiblast', 'hypoblast')] ->.;
    ## prettify stage
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, before.or.after.8.cell.ordered:=factor(before.or.after.8.cell, levels=unique(temp.stage.dt[, before.or.after.8.cell]))] ->.;
    ## 
    merge(x=.[, list(stage, Gene_ID, Gene_Name)], y=RE.valid.genes.only.dt[, list(CHROM, POS, stage, Gene_ID, Gene_Name, SAMPLE, gse, Annotation)],
          by=c("stage", "Gene_ID", "Gene_Name"),
          all.x=TRUE, all.y=FALSE)[, list(CHROM, POS, Gene_ID, Gene_Name, Annotation, SAMPLE, gse, stage)] -> .;
}
fwrite(postimplantation.RE.genes.and.their.RE.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.gene/postimplantation.RE.genes.and.their.RE.dt.csv.gz")

