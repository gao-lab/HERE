library("data.table")
library("magrittr")

merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename"]])

not.used.variable <- '
merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- fread("result/S51_5__filter_for_A_to_G_sites/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")
'

cat("dim of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt", dim(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt), "\n")

phenotype.output.at.gsm.level.dt <- fread(snakemake@input[["phenotype_output_at_gsm_level_dt_filename"]])
not.used.variable <- '
phenotype.output.at.gsm.level.dt <- fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt")
'


combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.zero.depth.records.only.dt <- fread(snakemake@input[["combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_zero_depth_records_only_txt_gz_filename"]]) %>% setnames(c("SAMPLE", "CHROM", "POS", "depth"))

not.used.variable <- '
## combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.zero.depth.records.only.dt <- fread("result/S52_2__concatenate_all_variant_coverages_of_merged_bam/201218-fifth-dataset/201221-fifth-phenotype-collection/combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.zero.depth.records.only.txt.gz")  %>% setnames(c("SAMPLE", "CHROM", "POS", "depth"))
'

cat("dim of combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.zero.depth.records.only.dt", dim(combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.zero.depth.records.only.dt), "\n")

test.merge.result.dt <- merge(x=merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt, y=combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.zero.depth.records.only.dt, by=c("SAMPLE", "CHROM", "POS"), all=FALSE)

cat("dim of test.merge.result.dt", dim(test.merge.result.dt), "\n")


fwrite(test.merge.result.dt, snakemake@output[["test_merge_result_dt_txt_filename"]])

rm(combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.zero.depth.records.only.dt)
gc()


combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.nonzero.depth.records.only.dt <- fread(snakemake@input[["combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_nonzero_depth_records_only_txt_gz_filename"]]) %>% setnames(c("SAMPLE", "CHROM", "POS", "depth"))

not.used.variable <- '
combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.nonzero.depth.records.only.dt <- fread("result/S52_2__concatenate_all_variant_coverages_of_merged_bam/201218-fifth-dataset/201221-fifth-phenotype-collection/combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.nonzero.depth.records.only.txt.gz")  %>% setnames(c("SAMPLE", "CHROM", "POS", "depth"))
'

cat("dim of combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.nonzero.depth.records.only.dt", dim(combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.nonzero.depth.records.only.dt), "\n")

sample.sequenced.dt <- rbindlist(lapply(readLines(paste(sep="", "external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/", snakemake@wildcards[['DATASET_COLLECTION_NAME']])), function(temp.DATASET_RNA_EDITING_NAME){return(unique(fread(paste(sep="", "external/DATASET_RNA_EDITING_NAME_DIRECTORY/", temp.DATASET_RNA_EDITING_NAME))))}))

cat("dim of sample.sequenced.dt", dim(sample.sequenced.dt), "\n")


not.used.variable <- '
sample.sequenced.dt <- rbindlist(lapply(readLines(paste(sep="", "external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/", "201218-fifth-dataset")), function(temp.DATASET_RNA_EDITING_NAME){return(unique(fread(paste(sep="", "external/DATASET_RNA_EDITING_NAME_DIRECTORY/", temp.DATASET_RNA_EDITING_NAME))))}))
'


merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt <- merge(x=merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, list(SAMPLE, CHROM, POS, ID, REF, ALT, event.summary, SUBSET, AC, AN, AF, site.occurrence)], y=combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.nonzero.depth.records.only.dt, by=c("SAMPLE", "CHROM", "POS"), all=TRUE) %>% {merge(x=., y=phenotype.output.at.gsm.level.dt, by.x="SAMPLE", by.y="gsm", all.x=TRUE, all.y=FALSE)} %>% {merge(x=., y=phenotype.output.at.gsm.level.dt[gsm %in% sample.sequenced.dt[, SAMPLE_NAME], list(total.sample.count.for.this.sample.group=.N), list(stage, is.normal)], by=c("stage", "is.normal"), all.x=TRUE, all.y=TRUE)}

cat("dim of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt", dim(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt), "\n")


fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt, snakemake@output[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_merged_with_coverage_dt_txt_gz_filename"]])

not.used.variable <- '
fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt, "result/S52_2__concatenate_all_variant_coverages_of_merged_bam/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt.txt.gz")
'
