library("data.table")
library("readxl")
library("glue")

GSE36552.gsms.vector <- c("GSM922146", "GSM922147", "GSM922148", "GSM922149", "GSM922150", "GSM922151", "GSM922152", "GSM922153", "GSM922154", "GSM922155", "GSM922156", "GSM922157", "GSM922158", "GSM922159", "GSM922160", "GSM922161", "GSM922162", "GSM922163", "GSM922164", "GSM922165", "GSM922166", "GSM922167", "GSM922168", "GSM922169", "GSM922170", "GSM922171", "GSM922172", "GSM922173", "GSM922174", "GSM922175", "GSM922176", "GSM922177", "GSM922178", "GSM922179", "GSM922180", "GSM922181", "GSM922182", "GSM922183", "GSM922184", "GSM922185", "GSM922186", "GSM922187", "GSM922188", "GSM922189", "GSM922190", "GSM922191", "GSM922192", "GSM922193", "GSM896803", "GSM896804", "GSM896805", "GSM896806", "GSM896807", "GSM896808", "GSM896809", "GSM896810", "GSM896811", "GSM896812", "GSM896813", "GSM896814")

all.GSE36552.results.dt <- rbindlist(lapply(GSE36552.gsms.vector, function(temp.gsm){
    fread(glue("result/PS81_2__Qiu2016_filter_variants/single-100/201104-GSE36552-full124-100/{temp.gsm}/__merged__/trim-15bp-off-3prime-and-filter-by-soapnuke/hg19.fa/Ensembl75.GRCh37/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-3.6.0/GATK-3.6.0-default/human.138.GRCh37/default/dbSNP138.and.1000GenomesPhase3.and.GoNL/RNA_editing.sites.annotation.splicing.filter.txt"))[, `:=`(GSE="GSE36552", GSM=temp.gsm)]
}), use.names=TRUE)

all.GSE36552.results.unique.position.bed.dt <- unique(all.GSE36552.results.dt[, list(Chromosome, start=Position-1, end=Position-1+1)])

all.GSE36552.results.unique.position.bed.filename <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/qiu2016.reimplemented.GSE36552.unique.position.bed"
fwrite(all.GSE36552.results.unique.position.bed.dt, all.GSE36552.results.unique.position.bed.filename, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

liftOver.directory <- glue("{all.GSE36552.results.unique.position.bed.filename}.liftOver.directory/")
system(glue("mkdir -p {liftOver.directory}"))

system(glue("for i in `seq 1 {nrow(all.GSE36552.results.unique.position.bed.dt)}`; do echo $i; tail -n +$i {all.GSE36552.results.unique.position.bed.filename} | head -n +1 > {liftOver.directory}/$i.bed; ./tools/liftOver {liftOver.directory}/$i.bed ./tools/hg19ToHg38.over.chain {liftOver.directory}/$i.liftOver.bed {liftOver.directory}/$i.unmapped.bed; done"))
## TODO waiting for its finishing



GSE44183.gsms.vector <- c("GSM1160112", "GSM1160113", "GSM1160114", "GSM1160115", "GSM1160116", "GSM1160117", "GSM1160118", "GSM1160119", "GSM1160120", "GSM1160121", "GSM1160122", "GSM1160123", "GSM1160124", "GSM1160125", "GSM1160126", "GSM1160127", "GSM1160128", "GSM1160129", "GSM1160138", "GSM1160139", "GSM1160140")


all.GSE44183.results.dt <- rbindlist(lapply(GSE44183.gsms.vector, function(temp.gsm){
    fread(glue("result/PS81_2__Qiu2016_filter_variants/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/{temp.gsm}/__merged__/filter-by-soapnuke/hg19.fa/Ensembl75.GRCh37/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-3.6.0/GATK-3.6.0-default/human.138.GRCh37/default/dbSNP138.and.1000GenomesPhase3.and.GoNL/RNA_editing.sites.annotation.splicing.filter.txt"))[, `:=`(GSE="GSE44183", GSM=temp.gsm)]
}), use.names=TRUE)

all.GSE44183.results.unique.position.bed.dt <- unique(all.GSE44183.results.dt[, list(Chromosome, start=Position-1, end=Position-1+1)])

all.GSE44183.results.unique.position.bed.filename <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/qiu2016.reimplemented.GSE44183.unique.position.bed"
fwrite(all.GSE44183.results.unique.position.bed.dt, all.GSE44183.results.unique.position.bed.filename, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

liftOver.directory <- glue("{all.GSE44183.results.unique.position.bed.filename}.liftOver.directory/")
system(glue("mkdir -p {liftOver.directory}"))

system(glue("for i in `seq 1 {nrow(all.GSE44183.results.unique.position.bed.dt)}`; do echo $i; tail -n +$i {all.GSE44183.results.unique.position.bed.filename} | head -n +1 > {liftOver.directory}/$i.bed; ./tools/liftOver {liftOver.directory}/$i.bed ./tools/hg19ToHg38.over.chain {liftOver.directory}/$i.liftOver.bed {liftOver.directory}/$i.unmapped.bed; done"))


all.GSE44183.results.unique.position.bed.liftOver.dt <- all.GSE44183.results.unique.position.bed.dt[, data.table(.SD, index=.I)][, {if (index %% 1000 == 0) {print(glue("{date()}: processing index {index}"))}; fread(glue("{liftOver.directory}/{index}.liftOver.bed"), header=FALSE, col.names=c("Chr.hg38", "start.hg38", "end.hg38"))}, list(index, Chr.hg19=Chromosome, start.hg19=start, end.hg19=end)]
all.GSE44183.results.unique.position.bed.liftOver.dt[, hg19.pos:=start.hg19+1]
all.GSE44183.results.unique.position.bed.liftOver.dt[, hg38.pos:=start.hg38+1]

fwrite(all.GSE44183.results.unique.position.bed.liftOver.dt, glue("{all.GSE44183.results.unique.position.bed.filename}.liftOver.bed"))

all.GSE44183.results.with.liftOver.dt <- merge(x=all.GSE44183.results.dt, y=all.GSE44183.results.unique.position.bed.liftOver.dt[, list(Chr.hg19, hg19.pos, Chr.hg38, hg38.pos)],
      by.x=c("Chromosome", "Position"),
      by.y=c("Chr.hg19", "hg19.pos"),
      all.x=TRUE, all.y=TRUE)


## per-sample distribution, too high compared to reported
## > unique(all.GSE44183.results.with.liftOver.dt[, list(Chr.hg38, hg38.pos, GSM)])[, .N, GSM][, summary(N)]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    1746    4543    8083    7326    8838   13285 

## diagnosis: check all p-value-based filters

all.GSE44183.MutDet.results.dt <- rbindlist(lapply(GSE44183.gsms.vector, function(temp.gsm){
    fread(glue("result/PS81_1__Qiu2016_call_variants/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/{temp.gsm}/__merged__/filter-by-soapnuke/hg19.fa/Ensembl75.GRCh37/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-3.6.0/GATK-3.6.0-default/human.138.GRCh37/default/variation.sites.txt.gz"))[, `:=`(GSE="GSE44183", GSM=temp.gsm)]
}), use.names=TRUE)
setnames(all.GSE44183.MutDet.results.dt, c("Chr", "Pos", "Depth", "Ref", "Ref.Sup", "Ref.Mid", "Ref.End", "Ref.Plus", "Ref.Minus", "Ref.Mism", "Ref.NoMism", "Alt", "Alt.Sup", "Alt.Mid", "Alt.End", "Alt.Plus", "Alt.Minus", "Alt.Mism", "Alt.NoMism", "Index", "rposPvalue", "sbPvalue", "mmPvalue", "GSE", "GSM"))

all.GSE44183.MutDet.results.with.Step123.info.dt <- {
    all.GSE44183.MutDet.results.dt -> .;
    ## Step 1
    .[, `:=`(Step.1.pvalue=mmPvalue, mapped.reads=Depth, variant.supporting.reads=Alt.Sup, mismatch.frequency=Alt.Sup/Depth)] -> .;
    ## Step 2
    .[, `:=`(Step.2.pvalue=sbPvalue, variant.strand.frequency=Alt.Plus/Alt.Sup, reference.strand.frequency=Ref.Plus/Ref.Sup)] -> .;
    .[, `:=`(variant.strand.preference=abs(variant.strand.frequency - 0.5), reference.strand.preference=abs(reference.strand.frequency - 0.5))] -> .;
    ## Step 3
    .[, `:=`(Step.3.pvalue=rposPvalue, read.end.frequency=Alt.End/Alt.Sup, read.middle.frequency=Alt.Mid/Alt.Sup)] -> .;
    .
}


all.GSE44183.results.with.liftOver.and.Step123.info.dt <- merge(x=all.GSE44183.results.with.liftOver.dt, y=all.GSE44183.MutDet.results.with.Step123.info.dt,
      by.x=c("Chromosome", "Position", "Reference", "Alteration", "GSE", "GSM"),
      by.y=c("Chr", "Pos", "Ref", "Alt", "GSE", "GSM"),
      all.x=TRUE, all.y=FALSE)

our.edits.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")
our.edits.GSE44183.dt <- our.edits.dt[gse %in% c("GSE44183")]

all.GSE44183.results.with.liftOver.and.Step123.info.and.overlap.with.ours.dt <- {
    copy(all.GSE44183.results.with.liftOver.and.Step123.info.dt) -> .;
    .[, detected.in.Qiu2016:=TRUE] -> .;
    .[, sample.occurrence:=length(unique(GSM)), list(Chromosome, Position, Reference, Alteration)] -> .;
    .[,`:=`(
        Step1.pvalue.keep=mmPvalue<0.01,
        Step1.mapped.reads.keep=mapped.reads>=4,
        Step1.variant.supporting.reads.keep=variant.supporting.reads>=3,
        Step1.mismatch.frequency.keep=mismatch.frequency>=0.1,
        Step2.pvalue.keep=((sbPvalue<0.005) & (variant.strand.frequency > reference.strand.frequency)) == FALSE,
        Step2.variant.strand.frequency.keep=((variant.strand.frequency > 0.9) | (variant.strand.frequency < 0.1)) == FALSE,
        Step3.pvalue.keep=((rposPvalue<0.05) & (read.end.frequency > read.middle.frequency)) == FALSE,
        Step3.read.end.frequency.keep=(read.end.frequency > 0.9) == FALSE,
        Step7.sample.occurrence.keep=sample.occurrence>=2,
        Step8.mismatch.frequency.not.too.high.keep=(mismatch.frequency>0.95)==FALSE
    )] -> .;
    merge(x=., y=data.table(our.edits.GSE44183.dt, detected.in.ours=TRUE),
          by.x=c("GSE", "GSM", "Chr.hg38", "hg38.pos"),
          by.y=c("gse", "SAMPLE", "CHROM", "POS"),
          all.x=TRUE, all.y=TRUE) -> .;
    .
}

## examine GSM1160125
##  > .[GSM=="GSM1160125"][paste(sep="", Reference, "->", Alteration) %in% c("A->G", "T->C") | detected.in.ours, .N, list(detected.in.Qiu2016, detected.in.ours)]
##    detected.in.Qiu2016 detected.in.ours    N
## 1:                TRUE               NA 1267
## 2:                  NA             TRUE  563
## 3:                TRUE             TRUE  155

## 1. Apply sample occurrence filter
## .[GSM=="GSM1160125"][paste(sep="", Reference, "->", Alteration) %in% c("A->G", "T->C") | detected.in.ours][Step7.sample.occurrence.keep %in% c(TRUE, NA)][, .N, list( detected.in.Qiu2016, detected.in.ours)]
##    detected.in.Qiu2016 detected.in.ours   N
## 1:                  NA             TRUE 563
## 2:                TRUE               NA  76
## 3:                TRUE             TRUE  60

## 2. Apply mismatch-frequency-not-to-high filter
## .[GSM=="GSM1160125"][paste(sep="", Reference, "->", Alteration) %in% c("A->G", "T->C") | detected.in.ours][Step7.sample.occurrence.keep %in% c(TRUE, NA)][Step8.mismatch.frequency.not.too.high.keep %in% c(TRUE, NA)][, .N, list(detected.in.Qiu2016, detected.in.ours)]
##    detected.in.Qiu2016 detected.in.ours   N
## 1:                  NA             TRUE 563
## 2:                TRUE               NA  17
## 3:                TRUE             TRUE   5

## 3. diagnose each variant
## > nrow(.[GSM=="GSM1160125"][paste(sep="", Reference, "->", Alteration) %in% c("A->G", "T->C") | detected.in.ours][Step7.sample.occurrence.keep %in% c(TRUE, NA)][Step8.mismatch.frequency.not.too.high.keep %in% c(TRUE, NA)][is.na(detected.in.ours)])
## [1] 17



merge(x=all.GSE44183.results.dt, y=all.GSE44183.MutDet.results.dt[, list(Chr, Pos, Ref, Alt, rposPvalue, sbPvalue, mmPvalue, )]

qiu2016.dt <- data.table(read_excel("./external/papers/10.1186/s12864-016-3115-2/12864_2016_3115_MOESM2_ESM.xlsx", sheet="Table.S1", skip=1) )
## transform to hg38 coordinates before comparison

qiu2016.bed.dt <- qiu2016.dt[, list(Chr, start=Pos-1, end=Pos-1+1)]

qiu2016.bed.dt.filename <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/qiu2016.bed"
fwrite(qiu2016.bed.dt, qiu2016.bed.dt.filename, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

liftOver.directory <- glue("{qiu2016.bed.dt.filename}.liftOver.directory/")
system(glue("mkdir -p {liftOver.directory}"))

system(glue("for i in `seq 1 {nrow(qiu2016.bed.dt)}`; do echo $i; tail -n +$i {qiu2016.bed.dt.filename} | head -n +1 > {liftOver.directory}/$i.bed; ./tools/liftOver {liftOver.directory}/$i.bed ./tools/hg19ToHg38.over.chain {liftOver.directory}/$i.liftOver.bed {liftOver.directory}/$i.unmapped.bed; done"))

qiu2016.bed.liftOver.dt <- qiu2016.bed.dt[, data.table(.SD, index=.I)][, fread(glue("{liftOver.directory}/{index}.liftOver.bed"), header=FALSE, col.names=c("Chr.hg38", "start.hg38", "end.hg38")), list(index, Chr.hg19=Chr, start.hg19=start, end.hg18=end)]
qiu2016.bed.liftOver.dt[, pos:=start.hg38+1]

input.edits.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")

input.edits.dt[gse %in% c("GSE44183", "GSE36552")][stage %in% c("oocyte.mature", "oocyte.MII", "pronucleus", "zygote", "2-cell", "4-cell", "8-cell", "morula")][AF<0.95][, .N, ID]



'
> merge(x=qiu2016.bed.liftOver.dt[, list(Chr.hg38, pos, is.qiu2016=TRUE)], y=input.edits.dt[gse %in% c("GSE44183", "GSE36552")][stage %in% c("oocyte.mature", "oocyte.MII", "pronucleus", "zygote", "2-cell", "4-cell", "8-cell", "morula")][AF<0.95][, list(CHROM, POS, gse, SAMPLE, stage, is.ours=TRUE)], by.x=c("Chr.hg38", "pos"), by.y=c("CHROM", "POS"), all.x=TRUE, all.y=TRUE)[, list(Chr.hg38, pos, is.ours, is.qiu2016)] %>% unique %>% {.[, .N, list(is.qiu2016, is.ours)]}
   is.qiu2016 is.ours     N
1:         NA    TRUE 70980
2:       TRUE      NA  2920
3:       TRUE    TRUE  5893

'

RE.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.dt.txt.gz")

'
> merge(x=qiu2016.bed.liftOver.dt[, list(Chr.hg38, pos, is.qiu2016=TRUE)], y=RE.dt[gse %in% c("GSE44183", "GSE36552")][stage %in% c("oocyte.mature", "oocyte.MII", "pronucleus", "zygote", "2-cell", "4-cell", "8-cell", "morula")][AF<0.95 & (is.na(AF)==FALSE)][, list(CHROM, POS, gse, SAMPLE, stage, is.ours=TRUE)], by.x=c("Chr.hg38", "pos"), by.y=c("CHROM", "POS"), all.x=TRUE, all.y=TRUE)[, list(Chr.hg38, pos, is.ours, is.qiu2016)] %>% unique %>% {.[, .N, list(is.qiu2016, is.ours)]}
   is.qiu2016 is.ours    N
1:         NA    TRUE 5973
2:       TRUE      NA 5943
3:       TRUE    TRUE 2870

'
