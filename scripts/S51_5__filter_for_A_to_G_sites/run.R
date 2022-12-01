library("data.table")
library("magrittr")

merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_dt_txt_gz_filename"]])

not.used.variable <- '
merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt <- fread("result/S51_4__filter_for_variants_with_enough_sample_support/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt.txt.gz")
'

cat("dim of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt ", dim(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt), "\n")

combined.merged.variant.only.snpEff.event.summary.dt <- fread(snakemake@input[["combined_merged_variant_only_snpEff_event_summary_dt_gz_filename"]])

not.used.variable <- '
combined.merged.variant.only.snpEff.event.summary.dt <- fread("result/S18_1__combine_annotations/201218-fifth-dataset/__merged__/base-quality-no-smaller-than-25/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/snpEff/basic/10000000/combined.merged.variant.only.snpEff.event.summary.dt.txt.gz")

'

cat("dim of combined.merged.variant.only.snpEff.event.summary.dt", dim(combined.merged.variant.only.snpEff.event.summary.dt), "\n")


merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt <- merge(x=merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt, y=combined.merged.variant.only.snpEff.event.summary.dt, by.x="ID", by.y="ID", all.x=TRUE, all.y=FALSE)

cat("dim of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt", dim(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt), "\n")


fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt, snakemake@output[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_with_event_summary_dt_txt_gz_filename"]])

not.used.variable <- '
fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt, "result/S51_5__filter_for_A_to_G_sites/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz")
'


merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt[event.summary %in% c("A>G", "A>G;T>C")]

cat("dim of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt", dim(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt), "\n")

cat("result of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, list(SAMPLE, stage, is.normal)] %>% unique %>% {.[, .N, list(is.normal)]}", merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, list(SAMPLE, stage, is.normal)] %>% unique %>% {.[, .N, list(is.normal)]} %>% print %>% as.null, "\n")

cat("result of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, unique(ID)] %>% length", merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, unique(ID)] %>% length, "\n")

cat("result of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, list(ID, event.summary)] %>% unique %>% {.[, .N, event.summary]}", merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, list(ID, event.summary)] %>% unique %>% {.[, .N, event.summary]} %>% print %>% as.null, "\n")

fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt, snakemake@output[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename"]])

not.used.variable <- '
fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt, "result/S51_5__filter_for_A_to_G_sites/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")
'


merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.bed.dt <- unique(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, list(CHROM, POS-1, POS, ID, ".", ".")])

cat("dim of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.bed.dt", dim(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.bed.dt), "\n")

fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.bed.dt, snakemake@output[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_bed_filename"]], row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

not.used.variable <- '
fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.bed.dt, "result/S51_5__filter_for_A_to_G_sites/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
'
