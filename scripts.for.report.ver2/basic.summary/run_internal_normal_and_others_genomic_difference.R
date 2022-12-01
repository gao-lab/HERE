library("data.table")
library("readxl")
library("glue")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/genomic.distribution/"
dir.create(output.directory, recursive=TRUE)

input.edits.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")

merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt <- fread("result/S51_6__get_snpEff_annotation_subset_of_filtered_result/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt.txt.gz")


## valid sites, A>G
{
    merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt -> .;
    .[event == 'A>G'] -> .;
    .
} -> snpEff.annotation.A.to.G.events.only.dt 

## valid sites, gene level (with Annotation collapsed)
{
    ##
    snpEff.annotation.A.to.G.events.only.dt -> .;
    ## collapse at gene level
    .[, list(CHROM, POS, Annotation, Gene_Name, Gene_ID)] -> .;
    unique(.) -> .;
    ## collapse different annotations
    .[, Annotation.pasted:=paste(collapse=";", Annotation %>% sort %>% unique), list(CHROM, POS, Gene_Name, Gene_ID)] -> .;
    ##
    ## correct the Annotation
    .[, Annotation.corrected:="DUII_variant"] -> .;
    cat(date(), "looping over terms...\n");
    for (term in c("intron_variant", "non_coding_transcript_exon_variant", "synonymous_variant", "stop_retained_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", 'missense_variant', "splice_region_variant", "splice_acceptor_variant", "splice_donor_variant", "5_prime_UTR_premature_start_codon_gain_variant", "stop_lost", "start_lost")){
        cat(date(), "  updating term: ", term, "\n");
        .[grepl(term, Annotation.pasted), Annotation.corrected:=term] -> .;
    }
    . -> ..snpEff.annotation.A.to.G.events.only.Annotation.collapsed.dt
    fwrite(..snpEff.annotation.A.to.G.events.only.Annotation.collapsed.dt, glue("{output.directory}/snpEff.annotation.A.to.G.events.only.Annotation.collapsed.dt.gz"))
    ..snpEff.annotation.A.to.G.events.only.Annotation.collapsed.dt
} -> snpEff.annotation.A.to.G.events.only.Annotation.collapsed.dt

## valid sites, site level (collapsed across genes)
{
    snpEff.annotation.A.to.G.events.only.Annotation.collapsed.dt -> .;
    unique(.[, list(CHROM, POS, Annotation)]) -> .;
    ## collapse across genes
    .[, Annotation.pasted.across.genes:=paste(collapse=";", sort(Annotation)), list(CHROM, POS)] -> .;
    
    ## correct the Annotation
    .[, Annotation.corrected.across.genes:="DUII_variant"] -> .;
    cat(date(), "looping over terms...\n");
    for (term in c("intron_variant", "non_coding_transcript_exon_variant", "synonymous_variant", "stop_retained_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", 'missense_variant', "splice_region_variant", "splice_acceptor_variant", "splice_donor_variant", "5_prime_UTR_premature_start_codon_gain_variant", "stop_lost", "start_lost")){
        cat(date(), "  updating term: ", term, "\n");
        .[grepl(term, Annotation.pasted.across.genes), Annotation.corrected.across.genes:=term] -> .;
    }
    unique(.[, list(CHROM, POS, Annotation.corrected.across.genes)]) -> .;
    . -> ..snpEff.annotation.A.to.G.events.only.Annotation.collapsed.across.genes.dt
    fwrite(..snpEff.annotation.A.to.G.events.only.Annotation.collapsed.across.genes.dt, glue("{output.directory}/snpEff.annotation.A.to.G.events.only.Annotation.collapsed.across.genes.dt.gz"))
    ..snpEff.annotation.A.to.G.events.only.Annotation.collapsed.across.genes.dt
} -> snpEff.annotation.A.to.G.events.only.Annotation.collapsed.across.genes.dt


## valid sites, site level (with Annotation collapsed across gene) and merged with per-sample stats



merge(
    x=input.edits.dt,
    y=snpEff.annotation.A.to.G.events.only.Annotation.collapsed.across.genes.dt,
    by=c("CHROM", "POS"), all.x=TRUE, all.y=FALSE
) -> input.edits.with.Annotation.collapsed.across.genes.dt

fwrite(input.edits.with.Annotation.collapsed.across.genes.dt, glue("{output.directory}/input.edits.with.Annotation.collapsed.across.genes.dt.gz"))


## per-stage count, barplot
{
    
    input.edits.with.Annotation.collapsed.across.genes.dt ->.
    ## simplify Annotation
    .[, Annotation.simplified:=c(
            "3_prime_UTR_variant"="3'-UTR",
            "5_prime_UTR_variant"="5'-UTR",
            "5_prime_UTR_premature_start_codon_gain_variant"="5'-UTR",
            "DUII_variant"="intergenic/intragenic",
            "intron_variant"="intron",
            "non_coding_transcript_exon_variant"= "exon (ncRNA)"
        )[Annotation.corrected.across.genes]] -> .;
    .[is.na(Annotation.simplified) == TRUE, Annotation.simplified:="splicing/coding sequence"] -> .;
    ## get sample group only
    unique(.[, list(ID, Annotation.simplified, stage, is.normal)]) -> .;
    ## keep stages that are comparable
    .[, `:=`(count.of.different.is.normal=length(unique(is.normal))), stage] -> .;
    .[count.of.different.is.normal == 2] -> .;
    ## redefine source
    .[, list(source.raw=paste(collapse=";", is.normal)), list(stage, ID, Annotation.simplified)] -> .;
    .[, source:=c("TRUE"="normal only", "FALSE"="abnormal only", "TRUE;FALSE"="both", "FALSE;TRUE"="both")[source.raw]] -> ..annotated.dt;
    
    ## get counts
    ..annotated.dt[, list(count=.N), list(Annotation.simplified, stage, source)] -> .;
    ## prettify stage
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ## prettify Annotation.simplified
    .[, Annotation.simplified.ordered:=factor(Annotation.simplified, levels=c("3'-UTR", "intron", "exon (ncRNA)", "5'-UTR", "splicing/coding sequence", "intergenic/intragenic"))] -> .;
    . -> ..to.plot.dt
    ##
    ## plot 
    ggplot(..to.plot.dt, aes(x=source, y=count, fill=Annotation.simplified.ordered)) ->.;
    . + geom_bar(stat="identity", position="fill") ->.;
    . + coord_flip() ->.;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "top") ->.;
    . + facet_grid(stage.description.ordered~.) -> .;
    . + guides(fill=guide_legend(nrow=3)) -> .;
    . + labs(x="", y="Ratio", fill="") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/genomic.distribution/genomic.difference.between.normal.and.abnormal.per.stage.png",
        plot=.,
        width.r=0.6, height.r=0.7
    )

    fwrite(..to.plot.dt, glue("{output.directory}/source.data.for.genomic.difference.between.normal.and.abnormal.per.stage.png.csv.gz"))
    
    ##
    melt(dcast(..to.plot.dt[, list(count=sum(count)), list(stage.description.ordered, source, Annotation.simplified.collapsed=c("intron"="3_prime_UTR_or_intron_variant", "3'-UTR"="3_prime_UTR_or_intron_variant", "intergenic/intragenic"="other_variant", "exon (ncRNA)"="other_variant", "splicing/coding sequence"="other_variant", "5'-UTR"="other_variant")[Annotation.simplified])], stage.description.ordered + Annotation.simplified.collapsed ~ source, value.var="count"), measure.vars=c("normal only", "abnormal only"), variable.name="group", value.name="group.count")[, list(chisq.test.pvalue.for.3UTRorintron.others.vs.both.unique=chisq.test(.SD[, list(both, group.count)])$p.value), list(group, stage.description.ordered)][, p.value.adjusted:=p.adjust(chisq.test.pvalue.for.3UTRorintron.others.vs.both.unique, method="BH")][, p.value.adjusted.clipped:=pmax(p.value.adjusted, 2.2e-16)] -> .;
    .[, group.prettified:=c("normal only"="both vs. normal only", "abnormal only"="both vs. abnormal only")[group]] -> .;
    ggplot(., aes(x=group.prettified, y=-1*log10(p.value.adjusted.clipped))) ->.;
    . + geom_bar(stat="identity") ->.;
    . + coord_flip() ->.;
    . + geom_hline(yintercept = -1*log10(0.1), color="red", linetype="dashed") -> .;
    . + theme_pubr() ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "top") ->.;
    . + facet_grid(stage.description.ordered~.) -> .;
    . + scale_y_continuous(limits=c(0, 16), breaks=c(0, 1, 5, 10, 15)) -> .;
    . + labs(x="", y="-1 * log10(BH-adjusted p-value of chi-square test on ratio of\n3'-UTR-or-intron edits, between the two groups indicated by the y-axis \n(upper-clipped by -log10(2.2e-16))", fill="") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/genomic.distribution/genomic.difference.between.normal.and.abnormal.per.stage.pvalue.bar.png",
        plot=.,
        width.r=0.9, height.r=0.7
    )
    ##
    
}

## additional: check sequence depth distribution between different sources
{

    input.edits.with.Annotation.collapsed.across.genes.dt -> .;
    ##
    .[, list(ID, stage, SAMPLE, AN)] -> .;
    ## we need to focus on only those stages where normal and abnormal samples both exists 
    merge(x=., y=..annotated.dt[, list(ID, stage, source)], by=c("ID", "stage"), all.x=FALSE, all.y=FALSE) -> .;
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    . -> ..to.plot.depth.dt

    ggplot(..to.plot.depth.dt, aes(x=source, y=AN)) ->.;
    . + geom_boxplot() ->.;
    ##. + stat_compare_means(comparisons=list(c("both", "normal only"), c("both", "abnormal only")), method="wilcox.test", method.args=list(alternative="greater")) -> .;
    . + coord_flip() ->.;
    . + theme_pubr() ->.;
    . + facet_grid(stage.description.ordered~.) -> .;
    . + scale_y_log10() -> .;
    . + labs(x="", y="sequencing depth (in log10 scale)", fill="") ->.;
    ggsave.A4(
        filename="./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/genomic.distribution/genomic.difference.between.normal.and.abnormal.per.stage.AN.boxplot.png",
        plot=.,
        width.r=0.9, height.r=0.7
    )
    
    
}
