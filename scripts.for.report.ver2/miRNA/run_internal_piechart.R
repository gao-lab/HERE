library("data.table")
library("readxl")
library("magrittr")
source("./scripts/common/ggpubr.A4.R")

observed.edits.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.observed.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

RE.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

edited.ts.human.compared.with.original.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/miRNA/edited.ts.human.compared.with.original.dt.csv.gz")

input.reference.GTF.filename <- "./external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf"
gene.transcript.mapping.dt <- fread(cmd=paste(sep="", "cat ", input.reference.GTF.filename, " | grep -v '^#' | grep -P '\ttranscript\t' | cut -f 1,9 | sed -E -e 's@^(chr[^\t]+)\t.*gene_id \"([^\"]+)\"; transcript_id \"([^\"]+)\";.*@\\1\t\\2\t\\3@' "), header=FALSE, sep="\t", col.names=c("CHROM", "gene.id", "transcript.id"))


edit.and.gene.and.MBS.affected.type.dt <- {
    edited.ts.human.compared.with.original.dt -> .;
    ##s restrict to on-miRNA-binding-site edits
    .[edit.rel.POS.wrt.3UTR >= UTR.start-1 & edit.rel.POS.wrt.3UTR <= UTR.end+1] -> .;
    ##m merge with gene-transcript
    merge(
        x=., y=gene.transcript.mapping.dt,
        by.x=c("transcript.id"), by.y="transcript.id", all.x=TRUE, all.y=FALSE
    ) ->.;
    ## For each [CHROM, edit.POS, gene.id, miRNA.site.affected.type], keep only the existence of miRNA and collapse to a single record
    ## `miRNA.site.affected.type.pasted` will be used for subtyping
    .[, list(miRNA.site.affected.type.pasted=miRNA.site.affected.type %>% sort %>% unique %>% paste(collapse=";")), list(CHROM, edit.POS, gene.id)] %>% unique -> .;
}


## get 3'-UTR subsets
generate.3UTR.subset.dt <- function(temp.dt){
    temp.dt -> .;
    ##m merge with edit.and.gene.and.MBS.affected.type.dt
    merge(
        x=edit.and.gene.and.MBS.affected.type.dt, y=.[Annotation=='3_prime_UTR_variant'],
        by.x=c("CHROM", "edit.POS", "gene.id"),
        by.y=c("CHROM", "POS", "Gene_ID"),
        all.x=TRUE, all.y=TRUE
    ) ->.;
    ## mark [transcript.id, edit.POS] without `miRNA.site.affected.type.pasted` as "not.miRNA.binding.site.related"
    .[is.na(miRNA.site.affected.type.pasted) == TRUE, miRNA.site.affected.type.pasted:="not.miRNA.binding.site.related"] ->.;
    ##= prettify
    .[, miRNA.site.affected.type.pasted.corrected:=c(
            "unchanged"="site unchanged",
            "gained"="site gained", "gained;unchanged"="site gained",
            "lost"="site lost", "lost;unchanged"="site lost",
            "gained;lost"="mixed effect", "gained;lost;unchanged"="mixed effect",
            "not.miRNA.binding.site.related"="no overlaps"
        )[miRNA.site.affected.type.pasted]]
}

observed.edits.valid.3UTR.only.dt <- generate.3UTR.subset.dt(observed.edits.valid.genes.only.dt)
RE.valid.3UTR.only.dt <- generate.3UTR.subset.dt(RE.valid.genes.only.dt)

fwrite(observed.edits.valid.genes.only.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/observed.edits.valid.genes.only.dt.csv.gz")
fwrite(RE.valid.genes.only.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/RE.valid.genes.only.dt.csv.gz")


plot.piechart.for.3UTR.edits <- function(temp.dt){
    copy(temp.dt) -> .;
    ## collapse to [CHROM, edit.POS, gene.id, miRNA.site.affected.type.pasted.corrected, stage], factoring out [SAMPLE]
    ## equivalent to accounting for all edits that are observed in >=1 sample
    .[, list(CHROM, edit.POS, gene.id, miRNA.site.affected.type.pasted.corrected, stage)] ->.; unique(.) ->.;
    ## compute the count
    .[, list(count=.N), list(miRNA.site.affected.type.pasted.corrected, stage)] ->.;
    ## focus on core stages only
    .[stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell")] ->.;
    ## prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## prettify Annotation
    .[, miRNA.site.affected.type.pasted.corrected.ordered:=factor(
            miRNA.site.affected.type.pasted.corrected,
            levels=rev(c("no overlaps", "site unchanged", "site gained", "site lost", "mixed effect")))] ->.;
    ggplot(., aes(x=1, y=count, fill=miRNA.site.affected.type.pasted.corrected.ordered)) +
        geom_bar(stat='identity', position='fill') +
        coord_polar(theta="y") +
        facet_wrap(~stage.description.ordered, nrow=1) ->.;
    . +
        theme_pubr(base_size=10) +
        theme(axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), legend.position="right") +
        labs(x="", y="", fill="Edit type") +
        scale_fill_manual(values=rev(c("grey50", "lightgreen", "#FF8270", "#7087FF", "#B42FFF"))) ->.;
    .
}

{
    observed.edit.valid.3UTR.only.ggplot <- plot.piechart.for.3UTR.edits(observed.edits.valid.3UTR.only.dt)
    RE.valid.3UTR.only.ggplot <- plot.piechart.for.3UTR.edits(RE.valid.3UTR.only.dt)
    combined.valid.3UTR.only.ggplot <- ggarrange(
        plotlist=list(
            observed.edit.valid.3UTR.only.ggplot + theme(legend.position="none", plot.margin = unit(c(0,0,0,0), "pt")) + ggtitle("All edits"),
            RE.valid.3UTR.only.ggplot + theme(legend.position="bottom", plot.margin = unit(c(0,0,0,0), "pt")) + labs(fill="") + guides(fill=guide_legend(nrow=2)) + ggtitle("REs") 
        ), nrow=2, align="v", heights=c(1, 1.55))
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/3UTR.types.piechart.png",
        plot=combined.valid.3UTR.only.ggplot,
        width.r=0.47, height.r=0.3
    )
}


