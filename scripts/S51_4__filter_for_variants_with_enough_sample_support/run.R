library("data.table")
library("magrittr")



merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_dt_txt_gz_filename"]])

not.used.variable <- '
merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt <- fread("result/S51_3__filter_for_variants_with_enough_read_support/201218-fifth-dataset/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt.txt.gz")
'

cat("dim of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt: ", dim(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt), "\n")


phenotype.output.at.gsm.level.dt <- fread(snakemake@input[["phenotype_output_at_gsm_level_dt_filename"]])
not.used.variable <- '
phenotype.output.at.gsm.level.dt <- fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt")
'
cat("dim of phenotype.output.at.gsm.level.dt: ", dim(phenotype.output.at.gsm.level.dt), "\n")



sample.sequenced.dt <- rbindlist(lapply(readLines(paste(sep="", "external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/", snakemake@wildcards[['DATASET_COLLECTION_NAME']])), function(temp.DATASET_RNA_EDITING_NAME){return(unique(fread(paste(sep="", "external/DATASET_RNA_EDITING_NAME_DIRECTORY/", temp.DATASET_RNA_EDITING_NAME))))}))

not.used.variable <- '
sample.sequenced.dt <- rbindlist(lapply(readLines(paste(sep="", "external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/", "201218-fifth-dataset")), function(temp.DATASET_RNA_EDITING_NAME){return(unique(fread(paste(sep="", "external/DATASET_RNA_EDITING_NAME_DIRECTORY/", temp.DATASET_RNA_EDITING_NAME))))}))
'

cat("dim of sample.sequenced.dt: ", dim(sample.sequenced.dt), "\n")


merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt <- merge(x=merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt, y=phenotype.output.at.gsm.level.dt[gsm %in% sample.sequenced.dt[, SAMPLE_NAME]], by.x="SAMPLE", by.y="gsm", all.x=TRUE, all.y=TRUE)[, `:=`(site.occurrence=.N), list(ID, stage, is.normal)]


cat("dim of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt ", dim(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt), "\n")


fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt, snakemake@output[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_dt_txt_gz_filename"]])
not.used.variable <- '
fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt, "result/S51_4__filter_for_variants_with_enough_sample_support/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt.txt.gz")
'


phenotype.of.samples.without.edits.before.sample.support.filter.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt[is.na(ID) == TRUE]
cat("dim of phenotype.of.samples.without.edits.before.sample.support.filter.dt ", dim(phenotype.of.samples.without.edits.before.sample.support.filter.dt), "\n")

fwrite(phenotype.of.samples.without.edits.before.sample.support.filter.dt, snakemake@output[["phenotype_of_samples_without_edits_before_sample_support_filter_dt_filename"]])

not.used.variable <- '
fwrite(phenotype.of.samples.without.edits.before.sample.support.filter.dt, "result/S51_4__filter_for_variants_with_enough_sample_support/201218-fifth-dataset/201221-fifth-phenotype-collection/phenotype.of.samples.without.edits.before.sample.support.filter.dt.txt")
'


merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt[(SUBSET == "Alu") | (SUBSET != "Alu" & site.occurrence >= 2)]

cat("dim of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt ", dim(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt), "\n")


fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt, snakemake@output[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_dt_txt_gz_filename"]])
not.used.variable <- '
fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt, "result/S51_4__filter_for_variants_with_enough_sample_support/201218-fifth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt.txt.gz")
'
