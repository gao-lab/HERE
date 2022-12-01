library("data.table")
library("readxl")
library("magrittr")
library("glue")
library("eulerr")

qiu2016.dt <- data.table(read_excel("./external/papers/10.1186/s12864-016-3115-2/12864_2016_3115_MOESM2_ESM.xlsx", sheet="Table.S1", skip=1) )
## transform to hg38 coordinates before comparison

qiu2016.bed.dt <- qiu2016.dt[, list(Chr, start=Pos-1, end=Pos-1+1)]

qiu2016.bed.dt.filename <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/qiu2016.bed"
fwrite(qiu2016.bed.dt, qiu2016.bed.dt.filename, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

liftOver.directory <- glue("{qiu2016.bed.dt.filename}.liftOver.directory/")
system(glue("mkdir -p {liftOver.directory}"))

system(glue("for i in `seq 1 {nrow(qiu2016.bed.dt)}`; do echo $i; tail -n +$i {qiu2016.bed.dt.filename} | head -n +1 > {liftOver.directory}/$i.bed; ./tools/liftOver {liftOver.directory}/$i.bed ./tools/hg19ToHg38.over.chain {liftOver.directory}/$i.liftOver.bed {liftOver.directory}/$i.unmapped.bed; done"))

qiu2016.bed.liftOver.dt <- qiu2016.bed.dt[, data.table(.SD, index=.I)][, fread(glue("{liftOver.directory}/{index}.liftOver.bed"), header=FALSE, col.names=c("Chr.hg38", "start.hg38", "end.hg38")), list(index, Chr.hg19=Chr, start.hg19=start, end.hg19=end)]
qiu2016.bed.liftOver.dt[, `:=`(Chr=Chr.hg38, pos=start.hg38+1)]

fwrite(qiu2016.bed.liftOver.dt, glue("{qiu2016.bed.dt.filename}.liftOver.directory/qiu2016.bed.lifOver.dt.gz"))

input.edits.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")


merged.dt <- merge(x=qiu2016.bed.liftOver.dt[, list(Chr.hg38, pos, is.qiu2016=TRUE)], y=input.edits.dt[gse %in% c("GSE44183", "GSE36552")][stage %in% c("oocyte.mature", "oocyte.MII", "pronucleus", "zygote", "2-cell", "4-cell", "8-cell", "morula")][AF<0.95][, list(CHROM, POS, gse, SAMPLE, stage, is.ours=TRUE)], by.x=c("Chr.hg38", "pos"), by.y=c("CHROM", "POS"), all.x=TRUE, all.y=TRUE)

fwrite(merged.dt, glue("{qiu2016.bed.dt.filename}.merged.with.ours.dt.gz"))

merged.dt[is.na(is.qiu2016) == TRUE, is.qiu2016:=FALSE]
merged.dt[is.na(is.ours) == TRUE, is.ours:=FALSE]

merged.edits.only.dt <- unique(merged.dt[, list(Chr.hg38, pos, `detected in\nQiu et al.`=is.qiu2016, `detected in\nours`=is.ours)])

png("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/direct.comparison.with.Qiu2016.png", width=12, height=8, units="cm", res=600)
print(plot(euler(merged.edits.only.dt[, list(`detected in\nQiu et al.`, `detected in\nours`)]), quantities=TRUE, legend=list(labels=c("detected in\nQiu et al.", "detected in\nours")), fills=c("#FF8F8F", "lightblue")))
dev.off()

'
> merge(x=qiu2016.bed.liftOver.dt[, list(Chr.hg38, pos, is.qiu2016=TRUE)], y=input.edits.dt[gse %in% c("GSE44183", "GSE36552")][stage %in% c("oocyte.mature", "oocyte.MII", "pronucleus", "zygote", "2-cell", "4-cell", "8-cell", "morula")][AF<0.95][, list(CHROM, POS, gse, SAMPLE, stage, is.ours=TRUE)], by.x=c("Chr.hg38", "pos"), by.y=c("CHROM", "POS"), all.x=TRUE, all.y=TRUE)[, list(Chr.hg38, pos, is.ours, is.qiu2016)] %>% unique %>% {.[, .N, list(is.qiu2016, is.ours)]}
   is.qiu2016 is.ours     N
1:         NA    TRUE 70980
2:       TRUE      NA  2920
3:       TRUE    TRUE  5893

'
