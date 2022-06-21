library("data.table")
library("glue")

output.directory <- "report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/sample.stat/"
dir.create(output.directory, recursive=TRUE)

embryo.samples.dt <- fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt")

fwrite(embryo.samples.dt[(gse=="GSE36552") | (srr.mean.avgspotlen >= 75*2)], glue("{output.directory}/all.2071.samples.csv"))
