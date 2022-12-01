library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("scales")
library("foreach")
library("iterators")
library("clusterProfiler")
library("org.Hs.eg.db")
source("./scripts/common/ggpubr.A4.R")

## read total edits in all samples
edit.info.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")

## RE and their genes
RE.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")


normal.RE.and.gene.per.stage.dt <- {
    RE.valid.genes.only.dt ->.;
    ##s keep matched stages only
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula')] ->.;
    ##= keep columns needed 
    .[, list(CHROM, POS, Gene_ID, Gene_Name, Annotation.corrected, stage)] %>% unique ->.;
}



edit.info.to.check.embryos.only.dt <- {
    edit.info.dt ->.;
    ##s select datasets of interest only
    .[gse %in% c("GSE95477", "GSE133854")] ->.;
    ##s select pathological/old-mother embryos only
    .[is.normal==FALSE | gse == "GSE95477"] -> .;
    ##s keep matched stages only
    .[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula')] ->.;
    ##= add columns for phenotype
    .[gse=="GSE95477", maternal.bin:=c("young", "old")[(as.integer(maternal.age) > 35)+1] ]
    .[gse=="GSE133854", disease.renamed:=c("androgenetic"="AG", "parthenogenetic"="PG")[disease]]
}

disease.or.old.mother.embryo.lost.RE.dt <- {
    copy(normal.RE.and.gene.per.stage.dt) -> .;
    ## check lost status in GSE95477, old-mother embryos
    .[, completely.lost.in.GSE95477.old.mother.embryos:=FALSE]
    .[
        stage == 'oocyte.GV' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE95477" & stage == "oocyte.GV" & maternal.bin=="old", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE95477.old.mother.embryos:=TRUE]
    .[
        stage == 'oocyte.MII' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE95477" & stage == "oocyte.MII" & maternal.bin=="old", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE95477.old.mother.embryos:=TRUE]
    ## check lost status in GSE133854, AG embryos
    .[, completely.lost.in.GSE133854.AG.embryos:=FALSE]
    .[
        stage == 'zygote' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE133854" & stage == "zygote" & disease.renamed == "AG", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE133854.AG.embryos:=TRUE]    
    .[
        stage == '2-cell' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE133854" & stage == "2-cell" & disease.renamed == "AG", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE133854.AG.embryos:=TRUE]    
    .[
        stage == '4-cell' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE133854" & stage == "4-cell" & disease.renamed == "AG", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE133854.AG.embryos:=TRUE]    
    .[
        stage == '8-cell' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE133854" & stage == "8-cell" & disease.renamed == "AG", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE133854.AG.embryos:=TRUE]    
    .[
        stage == 'morula' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE133854" & stage == "morula" & disease.renamed == "AG", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE133854.AG.embryos:=TRUE]    
    ## check lost status in GSE133854, PG embryos
    .[, completely.lost.in.GSE133854.PG.embryos:=FALSE]
    .[
        stage == 'zygote' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE133854" & stage == "zygote" & disease.renamed == "PG", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE133854.PG.embryos:=TRUE]    
    .[
        stage == '2-cell' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE133854" & stage == "2-cell" & disease.renamed == "PG", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE133854.PG.embryos:=TRUE]    
    .[
        stage == '4-cell' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE133854" & stage == "4-cell" & disease.renamed == "PG", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE133854.PG.embryos:=TRUE]    
    .[
        stage == '8-cell' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE133854" & stage == "8-cell" & disease.renamed == "PG", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE133854.PG.embryos:=TRUE]    
    .[
        stage == 'morula' &
        ((paste(sep="", CHROM, "@", POS)
            %in%
            edit.info.to.check.embryos.only.dt[gse=="GSE133854" & stage == "morula" & disease.renamed == "PG", paste(sep="", CHROM, "@", POS)])
            == FALSE),
        completely.lost.in.GSE133854.PG.embryos:=TRUE]
    ##s seleect REs lost in at least one disease/old-mother embryo
    .[completely.lost.in.GSE95477.old.mother.embryos | completely.lost.in.GSE133854.AG.embryos | completely.lost.in.GSE133854.PG.embryos] ->.;
}

fwrite(disease.or.old.mother.embryo.lost.RE.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/disease.or.old.mother.embryo.lost.RE.dt.csv.gz")


{
    
    copy(disease.or.old.mother.embryo.lost.RE.dt) -> .;
    ##= add CHROM.and.POS column
    .[, CHROM.and.POS:=paste(sep="", CHROM, "@", POS)]
    ## melt on lost type
    melt(data=.,
         id.vars=c("CHROM.and.POS", "Gene_ID", "Annotation.corrected", "stage"),
         measure.vars=c("completely.lost.in.GSE95477.old.mother.embryos", "completely.lost.in.GSE133854.AG.embryos", "completely.lost.in.GSE133854.PG.embryos"),
         variable.name="lost.type", value.name="really.lost") ->.;
    ##s select RE-lost records only
    .[really.lost == TRUE] -> .;
    ##c get count of REs and affected genes
    .[, list(count.of.RE=length(unique(CHROM.and.POS)), count.of.gene=length(unique(Gene_ID))), list(Annotation.corrected, stage, lost.type)] ->.;
    ## melt on count type   
    melt(data=., id.vars=c("Annotation.corrected", "stage", "lost.type"), variable.name="entity", value.name="count") ->.;
    ##s keep RE count only
    .[entity=="count.of.RE"] ->.;
    ##CJ fill missing combinations
    setkey(., Annotation.corrected, stage, lost.type)
    .[CJ(Annotation.corrected, stage, lost.type, unique=TRUE)] ->.;
    setnafill(., type="const", fill=0, cols="count")
    ## prettify stages
    temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.reverse.ordered:=factor(stage.description, levels=rev(temp.stage.dt[, stage.description]))] ->.;
    ##= prettify lost type
    .[, lost.type.renamed:=c("completely.lost.in.GSE95477.old.mother.embryos"="GSE95477\nembryos from\nold mothers", "completely.lost.in.GSE133854.AG.embryos"="GSE133854\nAG embryos", "completely.lost.in.GSE133854.PG.embryos"="GSE133854\nPG embryos")[lost.type]]
    ##= prettify Annotation
    .[, Annotation.corrected.renamed:=c("3_prime_UTR_variant"="3'-UTR", "intron_variant"="intronic", "missense_variant"="missense")[Annotation.corrected]]
    . -> to.plot.dt
    
    fwrite(to.plot.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/source.data.for.REE.disease.or.old.mother.embryo.lost.REE.count.barplot.png.csv.gz")
    
    ## start plotting
    to.plot.dt -> .;
    ggplot(., aes(x=stage.description.reverse.ordered, y=count, fill=Annotation.corrected.renamed)) ->.;
    . + geom_bar(stat="identity", position="dodge") ->.;
    . + coord_flip() ->.;
    ## facet
    . + facet_wrap(~lost.type.renamed, nrow=1) ->.;
    ## add theme
    . + theme_pubr() ->.;
    . + theme(legend.position="top") ->.;
    . + labs(x="", y="# REE-matching edits\ncompletely lost in this group", fill="") ->.;
    . + scale_y_continuous(breaks=c(0, 15, 30, 45)) -> .;
    ## save image
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/REE.disease.or.old.mother.embryo.lost.REE.count.barplot.png",
        plot=.,
        width.r=0.5, height.r=0.25)
}

{
    
    copy(disease.or.old.mother.embryo.lost.RE.dt) ->.;
    ## get no-gene-version Gene IDs
    .[, Gene_ID.no.version:=sub(pattern="\\.[0-9]+", replacement="", x=Gene_ID)]
    ## keep AG-lost edits only
    .[completely.lost.in.GSE133854.AG.embryos==TRUE] -> .;
    ## run GO analysis
    temp.all.enrichGO.results.list.collection <- {
        temp.subset.dt <- .
        temp.Gene_ID.no.version.results.vector <- unique(temp.subset.dt[, Gene_ID.no.version])
        temp.subset.enrichGO.results.list <- foreach(temp.ontology=c("BP", "MF", "CC")) %do% {
            enrichGO(temp.Gene_ID.no.version.results.vector, OrgDb='org.Hs.eg.db', keyType="ENSEMBL", ont=temp.ontology, pvalueCutoff=0.01)
        }
        list(temp.subset.enrichGO.results.list)
    }

    temp.all.enrichGO.results.combined.dt <- rbindlist(foreach(temp.subset.enrichGO.results.list=iter(temp.all.enrichGO.results.list.collection)) %do% {
        rbindlist(foreach(temp.ontology=c("BP", "MF", "CC"), temp.single.enrichGO.result=iter(temp.subset.enrichGO.results.list)) %do% {
            if (is.null(temp.single.enrichGO.result)){
                data.table(temp.single.combination.dt, ontology=temp.ontology)
            } else if (nrow(temp.single.enrichGO.result@result) == 0) {
                data.table(temp.single.combination.dt, ontology=temp.ontology)
            } else {
                data.table(data.table(temp.single.enrichGO.result@result), ontology=temp.ontology)
            }
        }, use.names=TRUE, fill=TRUE)
    }, use.names=TRUE, fill=TRUE)[is.na(ID) == FALSE][, p.adjusted.BH.used.in.the.manuscript:=p.adjust(p=pvalue, method="BH")]
    
}

fwrite(temp.all.enrichGO.results.combined.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/disease.and.old.mother.embryos.lost.RE.all.enrichGO.results.combined.dt.csv.gz")




{
    
    copy(temp.all.enrichGO.results.combined.dt) -> .;
    ##s select p.adjusted significant   
    .[p.adjusted.BH.used.in.the.manuscript<0.1] -> .;
    ##s select GO terms with large gene counts
    .[Count>=4] -> .;
    ##= prettify Description
    .[, Description.renamed:=Description]
    .[Description=="protein-containing complex localization", Description.renamed:="protein-containing\ncomplex localization"]
    .[Description=="establishment of RNA localization", Description.renamed:="establishment of\nRNA localization"]
    .[Description=="nucleobase-containing compound transport", Description.renamed:="nucleobase-containing\ncompound transport"]
    .[Description=="purine-containing compound metabolic process", Description.renamed:="purine-containing\ncompound metabolic process"]
    .[Description=="purine nucleotide metabolic process", Description.renamed:="purine nucleotide\nmetabolic process"]
    .[Description=="peptidyl-lysine modification", Description.renamed:="peptidyl-lysine\nmodification"]
    ## backup data
    temp.dt <- copy(.)
    
    fwrite(temp.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/source.data.for.REE.disease.or.old.mother.embryo.lost.REE.gene.GO.gene.count.ge.4.barplot.png.csv.gz")
    
    ## start plotting
    ggplot(., aes(x=Description.renamed, y=p.adjusted.BH.used.in.the.manuscript)) ->.;
    . + geom_bar(stat="identity") ->.;
    . + coord_flip() -> .;
    ## add theme
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(lineheight=0.75)) -> .;
    . + labs(x="Functions shared across\n>= 4 genes targeted by\nAG-lost REE-matching edits", y="BH-adjusted p-value   ") ->.;
    ## save image
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/REE.disease.or.old.mother.embryo.lost.REE.gene.GO.gene.count.ge.4.barplot.png",
        plot=.,
        width.r=0.45, height.r=0.25)    
}
