library("data.table")
library("glue")

output.directory <- "report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/AEI"
dir.create(output.directory, recursive=TRUE)


## determine per-step strand assignment

{

    {
        glue("{output.directory}/Alu.and.gencode.intersection") -> ..Alu.and.gencode.intersection.filename
        ## no strandness required 
        system(glue("bedtools intersect -wo -a external/UCSC.Table.Browser.repeatmasker/repFamily.Alu/repeatmasker.bed -b external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf > {..Alu.and.gencode.intersection.filename}"))
        
        ## read the intersection
        fread(
            ..Alu.and.gencode.intersection.filename, sep="\t", header=FALSE,
            col.names=c(
                c("Alu.chrom", "Alu.chromStart.zerobased.inclusive", "Alu.chromEnd.zerobased.exclusive", "Alu.name", "Alu.score", "Alu.strand"),
                c("Gene.annotation.seqname", "Gene.annotation.source", "Gene.annotation.feature",
                  "Gene.annotation.start.onebased.inclusive", "Gene.annotation.end.onebased.inclusive", "Gene.annotation.score", "Gene.annotation.strand", "Gene.annotation.frame",
                  "Gene.annotation.group"),
                "count.of.overlap.bases"
            )
        ) -> .;
        ##
        ## pick columns needed later
        .[, list(
            Alu.chrom=Alu.chrom,
            Alu.start=Alu.chromStart.zerobased.inclusive + 1,
            Alu.end=Alu.chromEnd.zerobased.exclusive + 1 - 1,
            Alu.name=Alu.name,
            Gene.annotation.feature,
            Gene.annotation.strand,
            Gene.annotation.group,
            count.of.overlap.bases
        )] -> .;
        .[, Alu.ID:=with(data=., expr=glue("{Alu.chrom}_{Alu.start}_{Alu.end}_{Alu.name}"))] -> .;
        
    } -> ..Alu.and.Gene.annotation.overlap.dt

    ## determine strand, step 1
    {
        copy(..Alu.and.Gene.annotation.overlap.dt) -> .;
        .[, list(step1.strand=paste(sep="", collapse=";", sort(unique(Gene.annotation.strand)))), list(Alu.ID)] -> .;
    } -> ..Alu.strand.annotation.step.1.dt
    fwrite(..Alu.strand.annotation.step.1.dt, glue("{output.directory}/Alu.strand.annotation.step.1.dt.gz"))
        
    ## determine strand, step 2
    {
        ## merge to get the step1-ambiguous records
        merge(x=..Alu.and.Gene.annotation.overlap.dt, y=..Alu.strand.annotation.step.1.dt,
              by="Alu.ID",
              all=FALSE) -> .;
        .[step1.strand == '-;+'] -> .;
        ## check overlaps with exons/transcripts
        .[, list(
            strands.overlapping.exons=paste(sep="", collapse=";", sort(unique(Gene.annotation.strand[Gene.annotation.feature == "exon"]))),
            strands.overlapping.transcripts=paste(sep="", collapse=";", sort(unique(Gene.annotation.strand[Gene.annotation.feature == "transcript"])))
        ), list(Alu.ID)] -> .;
        ## determine strand (step 2)
        .[strands.overlapping.exons == "-;+", step2.strand:="both.overlap.exons"] -> .;
        .[strands.overlapping.exons %in% c("+", "-"), step2.strand:=strands.overlapping.exons] -> .;
        .[(strands.overlapping.exons %in% c("-;+", "+", "-") == FALSE) & (strands.overlapping.transcripts %in% c("+", "-")), step2.strand:=strands.overlapping.transcripts] -> .;
        .[(strands.overlapping.exons %in% c("-;+", "+", "-") == FALSE) & (strands.overlapping.transcripts == "-;+"), step2.strand:="both.overlaps.introns"] -> .;
        ## add strand step 2
        .
    } -> ..Alu.strand.annotation.step.2.dt;
    fwrite(..Alu.strand.annotation.step.2.dt, glue("{output.directory}/Alu.strand.annotation.step.2.dt.gz"))

    ## [3] determine strand, step 3
    ## needs to account for per-sample expression
    
    ## [pre-3] get the expression table
    {        
        "result/BS06_1__get_expression_level/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/combined.gexpr.FPKM.matrix.txt" -> ..TEMP.FILENAME;
        fread(..TEMP.FILENAME, header=FALSE, skip=1) -> .
        setnames(., c("gene.id", scan(..TEMP.FILENAME, character(), nlines=1))) -> .;
        ## melt
        melt(data=., id.vars="gene.id", variable.name="SAMPLE", value.name="FPKM") -> .;
        ##
        .
    } -> ..combined.gexpr.FPKM.melt.dt
    
    ## [3-part1] get gene with max FPKM per sample
    {
        ##
        ## merge to get the step2-ambiguous records
        merge(x=..Alu.and.Gene.annotation.overlap.dt,
              y=..Alu.strand.annotation.step.2.dt,
              by="Alu.ID",
              all=FALSE) -> .;
        .[, Gene.annotation.gene.id:=sub(pattern='.*gene_id "([^"]+)".*', replacement="\\1", x=Gene.annotation.group)] -> .;
        unique(.[, list(Alu.ID, Gene.annotation.gene.id, Gene.annotation.strand)]) -> .;
        ## merge the Alu.ID ~ Gene id ~ Gene strand table with the per-sample expression table
        merge(
            x=., y=..combined.gexpr.FPKM.melt.dt,
            by.x="Gene.annotation.gene.id", by.y="gene.id",
            all.x=TRUE, all.y=FALSE,
            allow.cartesian=TRUE
        ) -> .;
        ## pick maximal FPKM record per SAMPLE x Alu.ID
        ## NOTE: if there are multiple maximal values, take all of them
        ##
        .[, max.FPKM:=max(FPKM), list(Alu.ID, SAMPLE)] -> .; ## ~3 minutes
        .[FPKM==max.FPKM] -> .;
        ##
        .
    } -> ..Alu.ID.and.gene.id.with.FPKM.dt
    fwrite(..Alu.ID.and.gene.id.with.FPKM.dt, glue("{output.directory}/Alu.ID.and.gene.id.with.FPKM.dt.gz"))

    ## [3 part-2] determine strand
    {
        ## The following strategy is way much faster than directly computing `paste(collapse=";")`
        unique(..Alu.ID.and.gene.id.with.FPKM.dt[, list(Alu.ID, SAMPLE, Gene.annotation.strand)]) -> .; ## ~20s
        setkey(., Alu.ID, SAMPLE) -> .; ## ~4s
        .[, list(step3.strand.count=.N, step3.strand.first=Gene.annotation.strand[1]), list(Alu.ID, SAMPLE)] -> .; ## ~10s
        .[step3.strand.count==1, step3.strand:=step3.strand.first] -> .;
        .[step3.strand.count==2, step3.strand:="-;+"] -> .;
        ##
        .
    } -> ..Alu.strand.annotation.step.3.dt
    fwrite(..Alu.strand.annotation.step.3.dt, glue("{output.directory}/Alu.strand.annotation.step.3.dt.gz"))

    ## step 4 needs the original AC-AN file
    ## will process it ad hoc
}


## Get the overlap between each Alu element and each edit
{
    ## [1] generate variant bed (A>G and C>T)
    {

        "result/S51_6__get_snpEff_annotation_subset_of_filtered_result/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.variant.only.snpEff.annotation.dt.txt.gz" -> .;
        fread(.) -> .;
        ##
        ## pick rows relevant to AEI and control index
        .[event %in% c("A>G", "A>G;T>C", "C>T", "C>T;G>A")] -> .;
        ## write annotation type for later validation
        fwrite(unique(.[, list(ID, Gene_ID, Feature_ID, Annotation)]), glue("{output.directory}/AEI.related.variants.snpEff.annotation.dt.gz"))
        ## focus on variants only
        unique(.[, list(ID, CHROM, POS)]) -> .;
        ## generate bed-format table
        .[, list(
            chrom=CHROM,
            chromStart=POS-1,
            chromEnd=POS-1+1,
            name=ID,
            score=0,
            strand="."           
        )] -> .;
        ##
        .
    } -> ..AEI.related.variants.bed.format.dt
    fwrite(..AEI.related.variants.bed.format.dt, glue("{output.directory}/AEI.related.variants.bed.format.dt.gz"))
    fwrite(..AEI.related.variants.bed.format.dt, glue("{output.directory}/AEI.related.variants.bed.format.bed"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    
    ## [2] get the Alu ~ variant table
    {
        ##
        glue("{output.directory}/Alu.bed.and.AEI.related.variants.bed.bedtools.intersect") -> .;        
        ## no strandness required in intersection
        system(glue("bedtools intersect -wa -wb -a external/UCSC.Table.Browser.repeatmasker/repFamily.Alu/repeatmasker.bed -b {output.directory}/AEI.related.variants.bed.format.bed > {.}")) ## shell command, return value is not used        
        ## read the intersection table
        fread(
            ., header=FALSE,
            col.names=c(
                c("Alu.chrom", "Alu.chromStart.bed", "Alu.chromEnd.bed", "Alu.name", "Alu.score", "Alu.strand"),
                c("variant.chrom", "variant.chromStart.bed", "variant.chromEnd.bed", "variant.ID", "variant.score", "variant.strand"))
        ) -> .;
        ## add Alu ID
        .[, Alu.ID:=with(data=., expr=glue("{Alu.chrom}_{Alu.chromStart.bed + 1}_{Alu.chromEnd.bed + 1 - 1}_{Alu.name}"))] -> .;
        ## pick relevant columns
        .[, list(
            Alu.ID,
            variant.ID        
        )] -> .;
        ##
        .
    } -> ..Alu.ID.and.AEI.related.variant.ID.overlap.dt
    fwrite(..Alu.ID.and.AEI.related.variant.ID.overlap.dt, glue("{output.directory}/Alu.ID.and.AEI.related.variant.ID.overlap.dt.gz"))
    ..Alu.ID.and.AEI.related.variant.ID.overlap.dt
    
} -> Alu.ID.and.AEI.related.variant.ID.overlap.dt



## Assign strands to each Alu per sample

{

    ## [1] build variant info per sample
    {
        ##
        "result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz" -> .;
        fread(.) -> .;
        ##
        ## pick relevant columns
        .[, list(
            variant.ID=ID,
            SAMPLE,
            AC,
            AN
        )] -> .;
        ##
        .
    } -> ..variant.ID.info.per.sample.dt

    ## [2] join with Alu - variant table
    {
        merge(x=..variant.ID.info.per.sample.dt, y=..Alu.ID.and.AEI.related.variant.ID.overlap.dt,
              by="variant.ID",
              all=FALSE) -> .;
        ##
        .
    } -> ..variant.ID.info.with.Alu.per.sample.dt

    ## [3] join strand step1 table
    {
        ##
        ..variant.ID.info.with.Alu.per.sample.dt -> .;
        ## merge strand 1 table
        merge(
            x=.,
            y=fread(glue("{output.directory}/Alu.strand.annotation.step.1.dt.gz")),
            by="Alu.ID",
            all.x=TRUE, all.y=FALSE
        ) -> .;
        ## the "NA" step1.strand values are due to their intergenic/upstream/downstream location (i.e., they are not within the gene body)
        ## merge(x=., y=fread(glue("{output.directory}/AEI.related.variants.snpEff.annotation.dt.gz")), by.x="variant.ID", by.y="ID", all.x=TRUE, all.y=FALSE, allow.cartesian = TRUE)[, .N, list(step1.strand, Annotation)][order(step1.strand)][is.na(step1.strand)]
        ##    step1.strand              Annotation      N
        ## 1:         <NA>       intergenic_region 514560
        ## 2:         <NA>   upstream_gene_variant 343413
        ## 3:         <NA> downstream_gene_variant 862587
        ## therefore, these hits should be discarded
        ##
        ## discard NA hits
        .[is.na(step1.strand) == FALSE] -> .;
        ##
        .
    } -> ..variant.ID.info.with.Alu.strand.step1.per.sample.dt;

    ## [4] join strand step2 table
    {
        
        ..variant.ID.info.with.Alu.strand.step1.per.sample.dt -> .;
        ## merge strand 2 table
        merge(
            x=.,
            y=fread(glue("{output.directory}/Alu.strand.annotation.step.2.dt.gz")),
            by="Alu.ID",
            all.x=TRUE, all.y=FALSE            
        ) -> .;
        ##
        .
    } -> ..variant.ID.info.with.Alu.strand.step1.step2.per.sample.dt

    ## [5] join strand step3 table
    {

        ..variant.ID.info.with.Alu.strand.step1.step2.per.sample.dt -> .;
        ## merge strand 3 table
        fread(glue("{output.directory}/Alu.strand.annotation.step.3.dt.gz")) -> ..temp.strand.step3.dt
        merge(
            x=.,
            y=..temp.strand.step3.dt,
            by=c("Alu.ID", "SAMPLE"),
            all.x=TRUE, all.y=FALSE
        ) -> .;        
    } -> ..variant.ID.info.with.Alu.strand.step1.step2.step3.per.sample.dt

    ## [6] compute strand step4 table
    {
        
        ..variant.ID.info.with.Alu.strand.step1.step2.step3.per.sample.dt -> .;
        ## pick those that cannot be determined in step 3
        .[step3.strand == "-;+"] -> .;
        ## mark group (AEI for A->G/T->C, or the control index for C->T/G->A) and genomic strand ('+' for A->G and C->T, and '-' for T->C and G->A)
        .[, genomic.plus.strand.variant:=sub(pattern=".*_([ACGT]{1}_[ACGT]{1}$)", replacement="\\1", x=variant.ID)] -> .;
        .[, group:=c("A_G"="AEI", "T_C"="AEI", "C_T"="control.index", "G_A"="control.index")[genomic.plus.strand.variant]] -> .;
        .[, genomic.strand:=c("A_G"="+", "T_C"="-", "C_T"="+", "G_A"="-")[genomic.plus.strand.variant]] -> .;
        
        ## compute sum(AC)/sum(AN) per Alu element
        dcast(., Alu.ID + SAMPLE + group ~ genomic.strand, value.var="variant.ID", fun.aggregate=length) -> .;
        .[`+` >= `-`, step4.strand := "+"] -> .;
        .[`+` < `-`, step4.strand := "-"] -> .;
        ##
        .
    } -> ..Alu.strand.annotation.step.4.dt
    fwrite(..Alu.strand.annotation.step.4.dt, glue("{output.directory}/Alu.strand.annotation.step.4.dt.gz"))

    ## [7] determine the final strand for each Alu element (x sample x group)
    {
        
        ..variant.ID.info.with.Alu.strand.step1.step2.step3.per.sample.dt -> .;
        ## get the group and genomic strand
        .[, genomic.plus.strand.variant:=sub(pattern=".*_([ACGT]{1}_[ACGT]{1}$)", replacement="\\1", x=variant.ID)] -> .;
        .[, group:=c("A_G"="AEI", "T_C"="AEI", "C_T"="control.index", "G_A"="control.index")[genomic.plus.strand.variant]] -> .;
        .[, genomic.strand:=c("A_G"="+", "T_C"="-", "C_T"="+", "G_A"="-")[genomic.plus.strand.variant]] -> .;
        ## join the step4 strand
        merge(
            x=.,
            y=..Alu.strand.annotation.step.4.dt[, list(Alu.ID, SAMPLE, group, step4.strand)],
            by=c("Alu.ID", "SAMPLE", "group"),
            all=TRUE
        ) -> .;
        ## determine the final strand
        ## (1) step1
        .[
        ((step1.strand %in% c("+", "-")) == TRUE),
        final.strand:=step1.strand
        ] -> .;
        ## (2) step2
        .[
        ((step1.strand %in% c("+", "-")) == FALSE) & ((step2.strand %in% c("+", "-")) == TRUE),
        final.strand:=step2.strand
        ] -> .; 
        ## (3) step3
        .[
        ((step1.strand %in% c("+", "-")) == FALSE) & ((step2.strand %in% c("+", "-")) == FALSE) & ((step3.strand %in% c("+", "-")) == TRUE),
        final.strand:=step3.strand
        ] -> .; 
        ## (4) step4
        .[
        ((step1.strand %in% c("+", "-")) == FALSE) & ((step2.strand %in% c("+", "-")) == FALSE) & ((step3.strand %in% c("+", "-")) == FALSE) & ((step4.strand %in% c("+", "-")) == TRUE),
        final.strand:=step4.strand
        ] -> .; 
        
        ## > .[, .N, list(group, step1.strand, step2.strand, step3.strand, step4.strand, final.strand)][order(paste(group, step1.strand, step2.strand, step3.strand, step4.strand, final.strand))]
        ##             group step1.strand          step2.strand step3.strand step4.strand final.strand       N
        ##  1:           AEI            +                  <NA>         <NA>         <NA>            + 1758500
        ##  2:           AEI            -                  <NA>         <NA>         <NA>            - 1764380
        ##  3:           AEI          -;+                     +            +         <NA>            +  132837
        ##  4:           AEI          -;+                     +            -         <NA>            +   34194
        ##  5:           AEI          -;+                     +          -;+            +            +     102
        ##  6:           AEI          -;+                     +          -;+            -            +      14
        ##  7:           AEI          -;+                     -            +         <NA>            -   22268
        ##  8:           AEI          -;+                     -            -         <NA>            -  123447
        ##  9:           AEI          -;+                     -          -;+            +            -      28
        ## 10:           AEI          -;+                     -          -;+            -            -      56
        ## 11:           AEI          -;+    both.overlap.exons            +         <NA>            +   55255
        ## 12:           AEI          -;+    both.overlap.exons            -         <NA>            -   49181
        ## 13:           AEI          -;+    both.overlap.exons          -;+            +            +      88
        ## 14:           AEI          -;+    both.overlap.exons          -;+            -            -      30
        ## 15:           AEI          -;+ both.overlaps.introns            +         <NA>            +   84762
        ## 16:           AEI          -;+ both.overlaps.introns            -         <NA>            -   86017
        ## 17:           AEI          -;+ both.overlaps.introns          -;+            +            +    2911
        ## 18:           AEI          -;+ both.overlaps.introns          -;+            -            -    3090
        ## 19: control.index            +                  <NA>         <NA>         <NA>            +   49035
        ## 20: control.index            -                  <NA>         <NA>         <NA>            -   47905
        ## 21: control.index          -;+                     +            +         <NA>            +    2956
        ## 22: control.index          -;+                     +            -         <NA>            +     889
        ## 23: control.index          -;+                     +          -;+            +            +       5
        ## 24: control.index          -;+                     +          -;+            -            +       5
        ## 25: control.index          -;+                     -            +         <NA>            -     783
        ## 26: control.index          -;+                     -            -         <NA>            -    2825
        ## 27: control.index          -;+                     -          -;+            +            -       4
        ## 28: control.index          -;+                     -          -;+            -            -       2
        ## 29: control.index          -;+    both.overlap.exons            +         <NA>            +    1203
        ## 30: control.index          -;+    both.overlap.exons            -         <NA>            -     948
        ## 31: control.index          -;+    both.overlap.exons          -;+            +            +       2
        ## 32: control.index          -;+    both.overlap.exons          -;+            -            -       2
        ## 33: control.index          -;+ both.overlaps.introns            +         <NA>            +    3948
        ## 34: control.index          -;+ both.overlaps.introns            -         <NA>            -    3922
        ## 35: control.index          -;+ both.overlaps.introns          -;+            +            +     226
        ## 36: control.index          -;+ both.overlaps.introns          -;+            -            -     198
        ##             group step1.strand          step2.strand step3.strand step4.strand final.strand       N

        ##
        .
    } -> ..Alu.strand.annotation.final.dt
    fwrite(..Alu.strand.annotation.final.dt, glue("{output.directory}/Alu.strand.annotation.final.dt.gz"))
    ## 
    ..Alu.strand.annotation.final.dt
} -> Alu.strand.annotation.final.dt


## Compute AEI

{
    ## [1] compute the AEI table
    {
        Alu.strand.annotation.final.dt -> .;
        ## only consider records with genomic.strand == final.strand
        .[genomic.strand == final.strand] -> .;
        ## compute the index
        .[, list(index.value=sum(AC)/sum(AN)*100), list(SAMPLE, group)] -> .;
        ## compute SNR
        dcast(., SAMPLE~group, value.var="index.value", fill=0)[, SNR:=AEI/control.index] -> .;
        ##
        .
    } -> ..AEI.dt
    fwrite(..AEI.dt, glue("{output.directory}/AEI.dt.gz"))
        
    ## [2] append ADAR expression
    ## [pre-2] get ADAR expression
    {
        ##
        "result/BS06_1__get_expression_level/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/combined.gexpr.FPKM.matrix.txt" -> ..TEMP.FILENAME;
        fread(..TEMP.FILENAME, header=FALSE, skip=1) -> .
        setnames(., c("gene.id", scan(..TEMP.FILENAME, character(), nlines=1))) -> .;
        ## ADARB1/ADAR2: ENSG00000197381.16
        ## ADAR1: ENSG00000160710.17
        .[gene.id %in% c("ENSG00000197381.16", "ENSG00000160710.17")] -> .;
        ## melt and dcast to get sample x gene
        melt(data=., id.vars="gene.id", variable.name="SAMPLE", value.name="FPKM") -> .;
        dcast(., SAMPLE ~ gene.id, value.var="FPKM") -> .;
        ##
        .        
    } -> ..ADAR.FPKM.dt
    fwrite(..ADAR.FPKM.dt, glue("{output.directory}/ADAR.FPKM.dt.gz"))

    ## [2-main] merge ADAR expression
    {
        merge(
            x=..AEI.dt, y=..ADAR.FPKM.dt,
            by="SAMPLE",
            all.x=TRUE, all.y=FALSE
        ) -> .;
        ##
        .
    } -> ..AEI.with.ADAR.FPKM.dt
    fwrite(..AEI.with.ADAR.FPKM.dt, glue("{output.directory}/AEI.with.ADAR.FPKM.dt.gz"))

    ## [3] append sample phenotype

    {
        fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt") -> .;
        merge(
            x=..AEI.with.ADAR.FPKM.dt, y=.,
            by.x="SAMPLE", by.y="gsm",
            all.x=TRUE, all.y=FALSE
        ) -> .;
        ##
        .
    } -> ..AEI.with.ADAR.FPKM.and.sample.info.dt
    fwrite(..AEI.with.ADAR.FPKM.and.sample.info.dt, glue("{output.directory}/AEI.with.ADAR.FPKM.and.sample.info.dt.gz"))

}
