library("data.table")
library("magrittr")

merged.variant.only.disjoint.with.population.variants.vcf.dt <- fread(cmd=paste(sep="", "zcat ", snakemake@input[["merged_variant_only_disjoint_with_population_variants_vcf_gz_filename"]], " | grep -v \"##\" "), header=TRUE) %>% setnames("#CHROM", "CHROM")

merged.long.dt <- fread(snakemake@input[["merged_long_table_dt_txt_gz_filename"]], header=TRUE)

not.used.variable <- '
merged.variant.only.disjoint.with.population.variants.vcf.dt <- fread(cmd=paste(sep="", "zcat ", "result/S51_2__filter_against_population_variants/201218-fifth-dataset/merged.variant.only.disjoint.with.population.variants.vcf.gz", "| grep -v \"##\" "), header=TRUE) %>% setnames("#CHROM", "CHROM")

merged.long.dt <- fread("result/S16_3__get_RNA_editing_site_long_table_from_a_dataset_collection/201218-fifth-dataset/__merged__/base-quality-no-smaller-than-25/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/merged.long.table.dt.txt.gz", header=TRUE)

'

cat("dim of merged.variant.only.disjoint.with.population.variants.vcf.dt: ", dim(merged.variant.only.disjoint.with.population.variants.vcf.dt), "\n")

cat("dim of merged.long.dt: ", dim(merged.long.dt), "\n")


merged.long.disjoint.with.population.dt <- merged.long.dt[ID %in% merged.variant.only.disjoint.with.population.variants.vcf.dt[, ID] == TRUE]

cat("dim of merged.long.disjoint.with.population.dt: ", dim(merged.long.disjoint.with.population.dt), "\n")

fwrite(merged.long.disjoint.with.population.dt, snakemake@output[["merged_long_disjoint_with_population_dt_txt_gz_filename"]])
not.used.variable <- '
fwrite(merged.long.disjoint.with.population.dt, "result/S51_3__filter_for_variants_with_enough_read_support/201218-fifth-dataset/merged.long.disjoint.with.population.dt.txt.gz")
'

merged.long.disjoint.with.population.ID.info.dt <- unique(merged.long.disjoint.with.population.dt[, list(ID, CHROM=sub(pattern="^([^_]+)_.*", replacement="\\1", x=ID), POS=sub(pattern="^[^_]+_([^_]+)_.*", replacement="\\1", x=ID))])[, count.of.types.of.mismatch:=.N, list(CHROM, POS)]

cat("dim of merged.long.disjoint.with.population.ID.info.dt: ", dim(merged.long.disjoint.with.population.ID.info.dt), "\n")

fwrite(merged.long.disjoint.with.population.ID.info.dt, snakemake@output[["merged_long_disjoint_with_population_dt_txt_gz_filename"]])
not.used.variable <- '
fwrite(merged.long.disjoint.with.population.ID.info.dt, "result/S51_3__filter_for_variants_with_enough_read_support/201218-fifth-dataset/merged.long.disjoint.with.population.ID.info.dt.txt.gz")
'


merged.long.disjoint.with.population.IDs.whose.coordinates.have.2.or.more.types.of.mismatch.vector <- merged.long.disjoint.with.population.ID.info.dt[count.of.types.of.mismatch >= 2, ID]

cat("length of merged.long.disjoint.with.population.IDs.whose.coordinates.have.2.or.more.types.of.mismatch.vector: ", length(merged.long.disjoint.with.population.IDs.whose.coordinates.have.2.or.more.types.of.mismatch.vector), "\n")


merged.long.disjoint.with.population.without.potential.polymorphism.dt <- merged.long.disjoint.with.population.dt[ (ID %in% merged.long.disjoint.with.population.IDs.whose.coordinates.have.2.or.more.types.of.mismatch.vector) == FALSE]

cat("dim of merged.long.disjoint.with.population.without.potential.polymorphism.dt: ", dim(merged.long.disjoint.with.population.without.potential.polymorphism.dt), "\n")


fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.dt, snakemake@output[["merged_long_disjoint_with_population_without_potential_polymorphism_dt_txt_gz_filename"]])
not.used.variable <- '
fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.dt, "result/S51_3__filter_for_variants_with_enough_read_support/201218-fifth-dataset/merged.long.disjoint.with.population.without.potential.polymorphism.dt.txt.gz")
'


merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.dt[(SUBSET=="Alu" & AN>=2 & AF>=0.1) | (SUBSET!="Alu" & AN>=2 & AF>=0.1 & AC>=2)]

cat("dim of merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt: ", dim(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt), "\n")

cat("length of unique(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt[, ID]): ", length(unique(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt[, ID])), "\n")

cat("length of unique(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt[, SAMPLE]): ", length(unique(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt[, SAMPLE])), "\n")



fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt, snakemake@output[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_dt_txt_gz_filename"]])
not.used.variable <- '
fwrite(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt, "result/S51_3__filter_for_variants_with_enough_read_support/201218-fifth-dataset/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt.txt.gz")
'
