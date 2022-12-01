library("data.table")
library("readxl")
library("magrittr")
library("glue")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA.intersection.all.edits/"
dir.create(output.directory, recursive=TRUE)


observed.edits.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.observed.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

RE.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

all.edited.intersection.of.TargetScan.and.miRanda.intersection.compared.with.original.annotated.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA.intersection.all.edits/all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.dt.gz")

input.reference.GTF.filename <- "./external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf"
gene.transcript.mapping.dt <- fread(cmd=paste(sep="", "cat ", input.reference.GTF.filename, " | grep -v '^#' | grep -P '\ttranscript\t' | cut -f 1,9 | sed -E -e 's@^(chr[^\t]+)\t.*gene_id \"([^\"]+)\"; transcript_id \"([^\"]+)\";.*@\\1\t\\2\t\\3@' "), header=FALSE, sep="\t", col.names=c("CHROM", "gene.id", "transcript.id"))


{
    copy(all.edited.intersection.of.TargetScan.and.miRanda.intersection.compared.with.original.annotated.dt) -> .;
    ##s restrict to on-miRNA-binding-site edits
    .[the.miRNA.family.in.the.union.of.original.and.edited.has.ge.1.site.overlapping.the.edit == TRUE] -> .;
    ##m merge with gene-transcript
    merge(
        x=., y=gene.transcript.mapping.dt,
        by.x=c("transcript.id"), by.y="transcript.id", all.x=TRUE, all.y=FALSE
    ) ->.;
    ## For each [CHROM, edit.POS, gene.id, miRNA.site.affected.type], keep only the existence of miRNA and collapse to a single record
    ## `miRNA.site.affected.type.pasted` will be used for subtyping
    .[, list(miRNA.site.affected.type.pasted=miRNA.site.affected.type %>% sort %>% unique %>% paste(collapse=";")), list(CHROM, edit.POS, gene.id)] %>% unique -> .;
} -> all.edit.and.gene.and.MBS.affected.type.dt


## get 3'-UTR subsets
generate.3UTR.subset.dt <- function(temp.dt){
    temp.dt -> .;
    ##m merge with edit.and.gene.and.MBS.affected.type.dt
    merge(
        x=all.edit.and.gene.and.MBS.affected.type.dt, y=.[Annotation=='3_prime_UTR_variant'],
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

fwrite(observed.edits.valid.genes.only.dt, glue("{output.directory}/observed.edits.valid.genes.only.dt.csv.gz"))
fwrite(RE.valid.genes.only.dt, glue("{output.directory}/RE.valid.genes.only.dt.csv.gz"))


plot.piechart.for.3UTR.edits <- function(temp.dt, plot.no.overlaps, return.to.plot.dt.only=FALSE){
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
    if (plot.no.overlaps == FALSE) {
        .[miRNA.site.affected.type.pasted.corrected != "no overlaps"] -> .;
    } else {
        ## do not discard no-overlaps
    }
    if ( return.to.plot.dt.only == TRUE) {
        return(.)
    }
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

## plot with no-overlaps
{
    observed.edit.valid.3UTR.only.ggplot <- plot.piechart.for.3UTR.edits(observed.edits.valid.3UTR.only.dt, plot.no.overlaps=TRUE)
    RE.valid.3UTR.only.ggplot <- plot.piechart.for.3UTR.edits(RE.valid.3UTR.only.dt, plot.no.overlaps=TRUE)
    combined.valid.3UTR.only.ggplot <- ggarrange(
        plotlist=list(
            observed.edit.valid.3UTR.only.ggplot + theme(legend.position="none", plot.margin = unit(c(0,0,0,0), "pt")) + ggtitle("All edits"),
            RE.valid.3UTR.only.ggplot + theme(legend.position="bottom", plot.margin = unit(c(0,0,0,0), "pt")) + labs(fill="") + guides(fill=guide_legend(nrow=2)) + ggtitle("REEs") 
        ), nrow=2, align="v", heights=c(1, 1.55))
    ggsave.A4(
        filename=glue("{output.directory}/REE.3UTR.types.piechart.with.no.overlaps.png"),
        plot=combined.valid.3UTR.only.ggplot,
        width.r=0.47, height.r=0.3
    )
}


## plot without no-overlaps
{
    observed.edit.valid.3UTR.only.ggplot <- plot.piechart.for.3UTR.edits(observed.edits.valid.3UTR.only.dt, plot.no.overlaps=FALSE)
    RE.valid.3UTR.only.ggplot <- plot.piechart.for.3UTR.edits(RE.valid.3UTR.only.dt, plot.no.overlaps=FALSE)
    ## compute 2x2 contigency table test
    observed.edit.valid.3UTR.only.to.plot.dt <- plot.piechart.for.3UTR.edits(observed.edits.valid.3UTR.only.dt, plot.no.overlaps=FALSE, return.to.plot.dt.only=TRUE)
    RE.valid.3UTR.only.to.plot.dt <- plot.piechart.for.3UTR.edits(RE.valid.3UTR.only.dt, plot.no.overlaps=FALSE, return.to.plot.dt.only=TRUE)

    {
        rbindlist(list(
            observed.edit.valid.3UTR.only.to.plot.dt[, list(source='observed.edit', total.count=sum(count)), list(stage.description.ordered, is.site.gained=(miRNA.site.affected.type.pasted.corrected=='site gained'))],
            RE.valid.3UTR.only.to.plot.dt[, list(source='RE', total.count=sum(count)),list(stage.description.ordered, is.site.gained=(miRNA.site.affected.type.pasted.corrected=='site gained'))]
        ), use.names=TRUE) -> .;
        . -> to.compute.p.value.dt
        to.compute.p.value.dt -> .;        
        dcast(., stage.description.ordered + is.site.gained ~ source, value.var="total.count")[, list(chisq.pvalue.raw=chisq.test(as.matrix(.SD[, list(RE, observed.edit)]))$p.value), list(stage.description.ordered)][, chisq.pvalue.adjusted:=p.adjust(chisq.pvalue.raw, method="BH")] -> .;
        .
    } -> pvalue.to.plot.dt


    merge(
        x=dcast(to.compute.p.value.dt, stage.description.ordered ~ paste(sep="", source, "__is.site.gained.IS.", is.site.gained, "__total.count"), value.var="total.count"),
        y=pvalue.to.plot.dt,
        by="stage.description.ordered",
        all.x=TRUE, all.y=TRUE         
    ) -> to.save.source.dt
    fwrite(to.save.source.dt, glue("{output.directory}/source.data.for.REE.3UTR.types.piechart.without.no.overlaps.png.csv.gz"))
    
    {
        ggplot(pvalue.to.plot.dt, aes(x="a", y="a")) -> .;
        . + geom_text(aes(label=format(chisq.pvalue.adjusted, digits=3, scientific=TRUE))) -> .;
        . + facet_wrap(~stage.description.ordered, nrow=1) -> .;
        . + theme_pubr(base_size=10) -> .;
        . + theme(axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank()) -> .;
        . + labs(x="", y="") -> .;
        .
    } -> pvalue.ggplot

    
    combined.valid.3UTR.only.ggplot <- ggarrange(
        plotlist=list(
            observed.edit.valid.3UTR.only.ggplot + theme(legend.position="none", plot.margin = unit(c(0,0,0,0), "pt")) + ggtitle("All MBS-overlapping edits"),
            RE.valid.3UTR.only.ggplot + theme(legend.position="bottom", plot.margin = unit(c(0,0,0,0), "pt")) + labs(fill="") + guides(fill=guide_legend(nrow=2)) + ggtitle("MBS-overlapping REEs"),
            pvalue.ggplot + ggtitle("P-value for chi-square test on frequency of MBS-\ngaining edits between MBS-overlapping REEs\nand all MBS-overlapping edits, BH-adjusted")
        ), nrow=3, align="v", heights=c(1, 1.55, 1))
    ggsave.A4(
        filename=glue("{output.directory}/REE.3UTR.types.piechart.without.no.overlaps.png"),
        plot=combined.valid.3UTR.only.ggplot,
        width.r=0.47, height.r=0.5
    )
}


