library("data.table")
library("magrittr")
source("./scripts/common/logger.R")


merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_with_event_summary_dt_txt_gz_filename"]])

not.used.variable <- '
merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt <- fread("result/S51_5__filter_for_A_to_G_sites/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz")
'


merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename"]])

not.used.variable <- '
merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- fread("result/S51_5__filter_for_A_to_G_sites/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")
'

combined.merged.vcf.reformatted.with.snpEff.ANN.split.annotation.dt <- fread(snakemake@input[["combined_merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_gz_filename"]])
'
combined.merged.vcf.reformatted.with.snpEff.ANN.split.annotation.dt <- fread("result/S18_1__combine_annotations/201218-fifth-dataset/__merged__/base-quality-no-smaller-than-25/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/snpEff/basic/10000000/combined.merged.variant.only.snpEff.ANN.single.match.dt.txt.gz")
' %>% as.null


setkey(combined.merged.vcf.reformatted.with.snpEff.ANN.split.annotation.dt, ID)


## all events
merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.variant.only.snpEff.annotation.dt <- combined.merged.vcf.reformatted.with.snpEff.ANN.split.annotation.dt[ID %in% merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt[, ID]]

report.expr(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.variant.only.snpEff.annotation.dt %>% dim)

report.expr(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.variant.only.snpEff.annotation.dt %>% {.[, ID]} %>% unique %>% length)


fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.variant.only.snpEff.annotation.dt, snakemake@output[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_with_event_summary_variant_only_snpEff_annotation_dt_txt_gz_filename"]])
'
fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.variant.only.snpEff.annotation.dt, "result/S51_6__get_snpEff_annotation_subset_of_filtered_result/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.variant.only.snpEff.annotation.dt.txt.gz")
' %>% as.null


## A-to-G_only

merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt <- combined.merged.vcf.reformatted.with.snpEff.ANN.split.annotation.dt[ID %in% merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, ID]]


report.expr(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt %>% dim)

report.expr(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt %>% {.[, ID]} %>% unique %>% length)


fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt, snakemake@output[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_variant_only_snpEff_annotation_dt_txt_gz_filename"]])
'
fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt, "result/S51_6__get_snpEff_annotation_subset_of_filtered_result/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt.txt.gz")
' %>% as.null



merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.with.snpEff.annotation.dt <- merge(x=merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt, y=merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt, by=c("CHROM", "POS", "REF", "ALT", "ID"), all=FALSE, allow.cartesian=TRUE)

report.expr(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.with.snpEff.annotation.dt %>% dim)

report.expr(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.with.snpEff.annotation.dt %>% {.[, ID]} %>% unique %>% length)




fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.with.snpEff.annotation.dt, snakemake@output[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_with_snpEff_annotation_dt_txt_gz_filename"]])
'
fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.with.snpEff.annotation.dt, "result/S51_6__get_snpEff_annotation_subset_of_filtered_result/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.with.snpEff.annotation.dt.txt.gz")
' %>% as.null
