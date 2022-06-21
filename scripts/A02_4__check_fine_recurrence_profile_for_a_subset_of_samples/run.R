library("data.table")
library("magrittr")
source("./scripts/common/logger.R")

merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_merged_with_coverage_dt_txt_gz_filename"]])

if (FALSE){
    merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt <- fread("result/S52_3__mark_unsequenced_editing_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt.txt.gz")
}


subset.name <- snakemake@wildcards[['SUBSET_NAME']]

if (FALSE){
    subset.name <- "all.normal.samples"
    subset.name <- "GSE133854.all"
}

report.expr(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt %>% dim)


## start subsetting

subset.dt <- NULL

if (subset.name == "all.normal.samples"){
    subset.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt[is.normal==TRUE]
    subset.dt[, group:=paste(sep="", stage, "@", is.normal)]
} else if (subset.name == "GSE133854.all"){
    subset.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt[gse=="GSE133854"]
    subset.dt[, group:=paste(sep="", stage, "@", disease)]
    subset.dt[disease %in% c("", NA) == TRUE, group:=paste(sep="", stage, "@", "biparental")]
} 

if (is.null(subset.dt) == TRUE){
    stop(paste0("Unsupported subset.name ", subset.name, "\n"))
}


## generate datasets needed by downstream analyses

## 1. subset.dt
fwrite(subset.dt, snakemake@output[["subset_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.dt.txt.gz"))
}

## 2. recurrence profile

subset.site.recurrence.comparison.dt <- subset.dt %>%
    {.[is.na(AC)==FALSE,
       list(site.occurrence.for.this.group=.N),
       list(CHROM, POS, group)]}

fwrite(subset.site.recurrence.comparison.dt, snakemake@output[["subset_site_recurrence_comparison_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.site.recurrence.comparison.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.site.recurrence.comparison.dt.txt.gz"))
}

## 3. recurrent edits (site + group only)




if (subset.name == 'all.normal.samples'){

    total.sample.count.for.normal.stages.dt <- fread(snakemake@input[["total_sample_count_for_normal_stages_dt_csv_filename"]])
    if (FALSE) {
        total.sample.count.for.normal.stages.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.sample.count.for.normal.stages.dt.csv")
    }

    total.sample.count.for.normal.groups.mapping.vector <- total.sample.count.for.normal.stages.dt %>% {temp.vector <- .[, total.sample.count]; names(temp.vector) <- .[, paste(sep="", stage, "@TRUE")]; temp.vector}


    subset.site.recurrence.comparison.dt[, total.sample.count:=total.sample.count.for.normal.groups.mapping.vector[group]]
    subset.site.recurrence.comparison.dt[, group.occurrence.pct:=site.occurrence.for.this.group / total.sample.count]
    
    subset.site.recurrence.comparison.recurrent.edits.only.dt <- subset.site.recurrence.comparison.dt[group.occurrence.pct >= 0.5]
}  else if (subset.name == 'GSE133854.all'){

    total.sample.count.for.GSE133854.all.dt <- fread(snakemake@input[["total_sample_count_for_GSE133854_all_dt_csv_filename"]])
    if (FALSE) {
        total.sample.count.for.GSE133854.all.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.sample.count.for.GSE133854.all.dt.csv")
    }

    total.sample.count.for.normal.groups.mapping.vector <- total.sample.count.for.GSE133854.all.dt %>% {temp.vector <- .[, total.sample.count]; names(temp.vector) <- .[, paste(sep="", stage, "@", disease)]; temp.vector}

    subset.site.recurrence.comparison.dt[, total.sample.count:=total.sample.count.for.normal.groups.mapping.vector[group]]

    subset.site.recurrence.comparison.dt[, group.occurrence.pct:=site.occurrence.for.this.group / total.sample.count]
    
    subset.site.recurrence.comparison.recurrent.edits.only.dt <- subset.site.recurrence.comparison.dt[group.occurrence.pct >= 0.5]
    
} 


fwrite(subset.site.recurrence.comparison.recurrent.edits.only.dt, snakemake@output[["subset_site_recurrence_comparison_recurrent_edits_only_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.site.recurrence.comparison.recurrent.edits.only.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.site.recurrence.comparison.recurrent.edits.only.dt.txt.gz"))
}

## 4. subset.dt with recurrent edits only 
## here the recurrency is group-specific (i.e., if an editing site is recurrent in group A but not in group B, then the final subset.dt will only have this site kept in group A samples but not group B samples)

subset.recurrent.edits.only.dt <- merge(x=subset.dt, y=subset.site.recurrence.comparison.recurrent.edits.only.dt, by=c("CHROM", "POS", "group"), all.x=FALSE)

fwrite(subset.recurrent.edits.only.dt, snakemake@output[["subset_recurrent_edits_only_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.recurrent.edits.only.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.recurrent.edits.only.dt.txt.gz"))
}

## 5. subset.dt with recurrence comparison info
## must use CHROM.and.POS rather than CHROM, POS; otherwise all combinations of CHROM and POS will appear

subset.site.recurrence.comparison.CJ.dt <- subset.site.recurrence.comparison.dt %>%
    ## mark recurrent edits
    {merge(x=., y=subset.site.recurrence.comparison.recurrent.edits.only.dt[, list(CHROM, POS, group, recurrence.type="recurrent")], by=c("CHROM", "POS", "group"), all=TRUE)} %>%
    ## mark `detected-but-not-recurrent` edits
    {.[is.na(recurrence.type) == TRUE, recurrence.type:="detected but not recurrent"]} %>%
    ## mark `not-detected` edits
    {.[, list(CHROM.and.POS=paste(sep="", CHROM, "_", POS), group, site.occurrence.for.this.group, total.sample.count, group.occurrence.pct, recurrence.type)]} %>%
    setkey(CHROM.and.POS, group) %>%
    {.[CJ(CHROM.and.POS, group, unique=TRUE)]} %>%
    {.[is.na(recurrence.type) == TRUE, recurrence.type:="not detected"]} %>%
    {.[, `:=`(
         CHROM=sub(pattern="_.*", replacement="", x=CHROM.and.POS),
         POS=sub(pattern=".*_", replacement="", x=CHROM.and.POS) %>% as.integer
     )]}

fwrite(subset.site.recurrence.comparison.CJ.dt, snakemake@output[["subset_site_recurrence_comparison_CJ_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.site.recurrence.comparison.CJ.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.site.recurrence.comparison.CJ.dt.txt.gz"))
}

## 6. subset.dt with recurrent edits only plus gene annotations

### 6.1. annotations for subsetted & valid sites (falls onto a gene, only A>G edits, only on protein-coding transcripts)

#### 6.1.1. all subsetted sites
merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_variant_only_snpEff_annotation_dt_txt_gz_filename"]])

if (FALSE) {
    merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt <- fread("result/S51_6__get_snpEff_annotation_subset_of_filtered_result/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt.txt.gz")
}

snpEff.annotation.for.subset.recurrent.edits.dt <- merge(
    x=subset.recurrent.edits.only.dt[, list(CHROM, POS)] %>% unique,
    y=merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt[, list(CHROM, POS, Annotation, Annotation_Impact, Gene_Name, Gene_ID, Feature_Type, Feature_ID, Transcript_BioType, Rank, HGVS.c, HGVS.p, `cDNA.pos / cDNA.length`, `CDS.pos / CDS.length`, `AA.pos / AA.length`, Distance, `ERRORS / WARNINGS / INFO`, event)],
    by=c("CHROM", "POS"), all.x=FALSE, all.y=FALSE)

fwrite(snpEff.annotation.for.subset.recurrent.edits.dt, snakemake@output[["snpEff_annotation_for_subset_recurrent_edits_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(snpEff.annotation.for.subset.recurrent.edits.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/snpEff.annotation.for.subset.recurrent.edits.dt.txt.gz"))
}

## 6.1.2. valid subsetted sites, transcript level
snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt <- snpEff.annotation.for.subset.recurrent.edits.dt %>%
    {.[Annotation %in% c("downstream_gene_variant", "intergenic_region", "upstream_gene_variant") == FALSE]} %>%
    {.[event == 'A>G']} %>%
    {.[Transcript_BioType=='protein_coding']}

fwrite(snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt, snakemake@output[["snpEff_annotation_for_subset_recurrent_edits_on_valid_transcripts_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt.txt.gz"))
}


## 6.1.3. valid subsetted sites, gene level (with Annotation collapsed)
snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt <- snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt %>%
    {.[, list(CHROM, POS, Annotation, Gene_Name, Gene_ID)]} %>%
    unique %>%
    {
        cat(date(), "generating Annotation.pasted ...\n");
        .[, Annotation.pasted:=paste(collapse=";", Annotation %>% sort %>% unique), list(CHROM, POS, Gene_Name, Gene_ID)]
    } %>%
    {.[, Annotation.corrected:="intron_variant"]} %>%
    {
        cat(date(), "looping over terms...\n");
        for (term in c("synonymous_variant", "stop_retained_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", 'missense_variant', "splice_region_variant", "splice_acceptor_variant", "splice_donor_variant", "5_prime_UTR_premature_start_codon_gain_variant", "stop_lost", "start_lost")){
            cat(date(), "  updating term: ", term, "\n");
            .[grepl(term, Annotation.pasted), Annotation.corrected:=term]
        }
        .
    } %>%
    {
        .[, Annotation.class:="exonic.or.splicing.related"]
        .[Annotation.corrected=='intron_variant', Annotation.class:="purely.intronic"]
    }

## each CHROM+POS+Gene_Name+Gene_ID has only one possible Annotation class
if (FALSE) {
    snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt[, list(CHROM, POS, Gene_Name, Gene_ID, Annotation.class)] %>% unique %>% {.[, .N, list(CHROM, POS, Gene_Name, Gene_ID)][N!=1]}
'
    Empty data.table (0 rows and 5 cols): CHROM,POS,Gene_Name,Gene_ID,N
'
}

fwrite(snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt, snakemake@output[["snpEff_annotation_for_subset_recurrent_edits_on_valid_genes_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt.txt.gz"))
}



## 6.1.4. valid subsetted sites, gene level (with Annotation collapsed) and merged with per-sample stats

subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- merge(
    x=subset.recurrent.edits.only.dt,
    y=snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt,
    by=c("CHROM", "POS"), all=FALSE)


if (FALSE) {
    subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt[is.na(AC)==FALSE, .N, list(SAMPLE, stage, Annotation.class)] %>% {dcast(., SAMPLE+stage~Annotation.class, value.var="N")[, log2.ES.I.ratio:=log2(exonic.or.splicing.related/purely.intronic)]} %>% {.[, quantile(log2.ES.I.ratio, na.rm=TRUE) %>% t %>% data.table, list(stage)]}
    '
                       stage         0%        25%       50%       75%      100%
 1:            oocyte.mature  1.7266558  1.7381526 1.7496495 1.7564444 1.7632393
 2:               pronucleus  1.6610344  1.7198921 1.7818276 1.8080704 1.8154629
 3:                   zygote  1.8327470  1.9631347 1.9988545 2.0321900 2.1854928
 4:                   2-cell  1.6411453  1.7659527 1.8083454 1.8536289 2.0549083
 5:                   4-cell  1.3356030  2.1901838 2.2571578 2.3476338 2.8050110
 6:                   8-cell -0.2223924  1.6780719 1.9068906 2.1154772 3.4594316
 7:                   morula  0.0000000  1.0000000 1.5849625 2.0000000 2.8073549
 8: post-implantation embryo  1.2925699  1.3492237 1.4058775 1.4408383 1.4757991
 9:                     hESC  1.5849625  3.0223678 3.0660892 3.4150375 4.7548875
10:               zygote.2PN  1.6508995  1.6984058 1.7319101 1.7528945 1.8109662
11:               oocyte.MII  1.8329694  1.9418434 1.9923449 2.0476858 2.3582603
12:         blastocyst.early -2.0000000 -0.7075187 0.5849625 1.1961587 1.8073549
13:        blastocyst.middle  1.3785116  1.6264904 1.8744691 2.1333933 2.3923174
14:          blastocyst.late         NA         NA        NA        NA        NA
15:                      ICM -1.0000000  0.0000000 0.0000000 0.0000000 1.5849625
16:                       TE  1.7004397  1.8930848 2.0506261 2.2801079 2.7224660
17:                oocyte.GV  1.4785732  1.6477370 1.7056032 1.7379000 1.7732359
18:                oocyte.MI  1.4687301  1.4689737 1.4692172 1.4694607 1.4697043
19:               2/4/8-cell  1.2181169  1.2433917 1.2686666 1.2939414 1.3192163
20:               2-cell_3PN  0.9139172  0.9139172 0.9139172 0.9139172 0.9139172
21:               4-cell_3PN  0.8981112  0.9255977 0.9530842 0.9805708 1.0080573
22:                      STB -2.0000000 -0.4150375 0.0000000 0.5849625 2.0000000
23:                      CTB  0.0000000  1.0000000 1.5849625 1.5849625 1.5849625
24:                      MTB  1.0000000  2.0000000 2.3571228 2.6990625 3.8073549
25:                 epiblast         NA         NA        NA        NA        NA
26:                hypoblast         NA         NA        NA        NA        NA
27:                      EVT  1.1699250  1.6715704 2.1039464 2.5077657 4.2479275
                       stage         0%        25%       50%       75%      100%

'
}


fwrite(subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt, snakemake@output[["subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz"))
}


## 6.2. all observed edits in the subset

#### 6.2.1. all subsetted sites
#### note that for 'subset.dt', some edits might be completely absent in all samples
#### therefore we need to remove these edits by ourselves
#### for 6.1.1, this has been done in the step generating "recurrence comparison" table

subset.observed.edits.only.dt <- subset.dt[is.na(AC)==FALSE]

snpEff.annotation.for.subset.observed.edits.dt <- merge(
    x=subset.observed.edits.only.dt[, list(CHROM, POS)] %>% unique,
    y=merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt[, list(CHROM, POS, Annotation, Annotation_Impact, Gene_Name, Gene_ID, Feature_Type, Feature_ID, Transcript_BioType, Rank, HGVS.c, HGVS.p, `cDNA.pos / cDNA.length`, `CDS.pos / CDS.length`, `AA.pos / AA.length`, Distance, `ERRORS / WARNINGS / INFO`, event)],
    by=c("CHROM", "POS"), all.x=FALSE, all.y=FALSE)

fwrite(snpEff.annotation.for.subset.observed.edits.dt, snakemake@output[["snpEff_annotation_for_subset_observed_edits_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(snpEff.annotation.for.subset.observed.edits.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/snpEff.annotation.for.subset.observed.edits.dt.txt.gz"))
}

## 6.2.2. valid subsetted sites, transcript level
snpEff.annotation.for.subset.observed.edits.on.valid.transcripts.dt <- snpEff.annotation.for.subset.observed.edits.dt %>%
    {.[Annotation %in% c("downstream_gene_variant", "intergenic_region", "upstream_gene_variant") == FALSE]} %>%
    {.[event == 'A>G']} %>%
    {.[Transcript_BioType=='protein_coding']}

fwrite(snpEff.annotation.for.subset.observed.edits.on.valid.transcripts.dt, snakemake@output[["snpEff_annotation_for_subset_observed_edits_on_valid_transcripts_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(snpEff.annotation.for.subset.observed.edits.on.valid.transcripts.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/snpEff.annotation.for.subset.observed.edits.on.valid.transcripts.dt.txt.gz"))
}


## 6.2.3. valid subsetted sites, gene level (with Annotation collapsed)
snpEff.annotation.for.subset.observed.edits.on.valid.genes.dt <- snpEff.annotation.for.subset.observed.edits.on.valid.transcripts.dt %>%
    {.[, list(CHROM, POS, Annotation, Gene_Name, Gene_ID)]} %>%
    unique %>%
    {
        cat(date(), "generating Annotation.pasted ...\n");
        .[, Annotation.pasted:=paste(collapse=";", Annotation %>% sort %>% unique), list(CHROM, POS, Gene_Name, Gene_ID)]
    } %>%
    {.[, Annotation.corrected:="intron_variant"]} %>%
    {
        cat(date(), "looping over terms...\n");
        for (term in c("synonymous_variant", "stop_retained_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", 'missense_variant', "splice_region_variant", "splice_acceptor_variant", "splice_donor_variant", "5_prime_UTR_premature_start_codon_gain_variant", "stop_lost", "start_lost")){
            cat(date(), "  updating term: ", term, "\n");
            .[grepl(term, Annotation.pasted), Annotation.corrected:=term]
        }
        .
    } %>%
    {
        .[, Annotation.class:="exonic.or.splicing.related"]
        .[Annotation.corrected=='intron_variant', Annotation.class:="purely.intronic"]
    }

## each CHROM+POS+Gene_Name+Gene_ID has only one possible Annotation class
if (FALSE) {
    snpEff.annotation.for.subset.observed.edits.on.valid.genes.dt[, list(CHROM, POS, Gene_Name, Gene_ID, Annotation.class)] %>% unique %>% {.[, .N, list(CHROM, POS, Gene_Name, Gene_ID)][N!=1]}
'
    Empty data.table (0 rows and 5 cols): CHROM,POS,Gene_Name,Gene_ID,N
'
}

fwrite(snpEff.annotation.for.subset.observed.edits.on.valid.genes.dt, snakemake@output[["snpEff_annotation_for_subset_observed_edits_on_valid_genes_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(snpEff.annotation.for.subset.observed.edits.on.valid.genes.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/snpEff.annotation.for.subset.observed.edits.on.valid.genes.dt.txt.gz"))
}



## 6.2.4. valid subsetted sites, gene level (with Annotation collapsed) and merged with per-sample stats

subset.observed.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- merge(
    x=subset.observed.edits.only.dt,
    y=snpEff.annotation.for.subset.observed.edits.on.valid.genes.dt,
    by=c("CHROM", "POS"), all=FALSE)




fwrite(subset.observed.edits.only.with.snpEff.annotation.on.valid.genes.only.dt, snakemake@output[["subset_observed_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.observed.edits.only.with.snpEff.annotation.on.valid.genes.only.dt, paste(sep="", "result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.observed.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz"))
}
