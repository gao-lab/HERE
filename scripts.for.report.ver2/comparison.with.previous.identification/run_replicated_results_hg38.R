library("data.table")
library("readxl")
library("glue")
library("doMC")
library("foreach")
library("iterators")
library("stringr")
library("ggpubr")
library("eulerr")
source("./scripts/common/ggpubr.A4.R")


## 1. Read the sample IDs

comparison.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/hg38.based.comparison.with.Qiu2016/"

dir.create(comparison.directory, recursive=TRUE)

GSE36552.gsms.vector <- c("GSM922146", "GSM922147", "GSM922148", "GSM922149", "GSM922150", "GSM922151", "GSM922152", "GSM922153", "GSM922154", "GSM922155", "GSM922156", "GSM922157", "GSM922158", "GSM922159", "GSM922160", "GSM922161", "GSM922162", "GSM922163", "GSM922164", "GSM922165", "GSM922166", "GSM922167", "GSM922168", "GSM922169", "GSM922170", "GSM922171", "GSM922172", "GSM922173", "GSM922174", "GSM922175", "GSM922176", "GSM922177", "GSM922178", "GSM922179", "GSM922180", "GSM922181", "GSM922182", "GSM922183", "GSM922184", "GSM922185", "GSM922186", "GSM922187", "GSM922188", "GSM922189", "GSM922190", "GSM922191", "GSM922192", "GSM922193", "GSM896803", "GSM896804", "GSM896805", "GSM896806", "GSM896807", "GSM896808", "GSM896809", "GSM896810", "GSM896811", "GSM896812", "GSM896813", "GSM896814")
GSE44183.gsms.vector <- c("GSM1160112", "GSM1160113", "GSM1160114", "GSM1160115", "GSM1160116", "GSM1160117", "GSM1160118", "GSM1160119", "GSM1160120", "GSM1160121", "GSM1160122", "GSM1160123", "GSM1160124", "GSM1160125", "GSM1160126", "GSM1160127", "GSM1160128", "GSM1160129", "GSM1160138", "GSM1160139", "GSM1160140")



## 2. Read # uniquely mapped bases and select valid samples

GSE.and.GSM.dt <- rbindlist(list(
    data.table(GSE="GSE36552", GSM=GSE36552.gsms.vector),
    data.table(GSE="GSE44183", GSM=GSE44183.gsms.vector)
), use.names=TRUE)


all.GSE36552.number.of.uniquely.mapped.bases.dt <- rbindlist(lapply(GSE36552.gsms.vector, function(temp.gsm){
    data.table(
        GSE="GSE36552",
        GSM=temp.gsm,
        number.of.uniquely.mapped.bases=scan(glue("result/PS31_1__get_number_of_uniquely_mapped_bases/single-100/201104-GSE36552-full124-100/{temp.gsm}/__merged__/trim-15bp-off-3prime-and-filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-2.8-1/GATK-2.8-1-default/151/common_all/number.of.uniquely.mapped.bases"), what=numeric()))
}), use.names=TRUE)

all.GSE44183.number.of.uniquely.mapped.bases.dt <- rbindlist(lapply(GSE44183.gsms.vector, function(temp.gsm){
    data.table(
        GSE="GSE44183",
        GSM=temp.gsm,
        number.of.uniquely.mapped.bases=scan(glue("result/PS31_1__get_number_of_uniquely_mapped_bases/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/{temp.gsm}/__merged__/filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-2.8-1/GATK-2.8-1-default/151/common_all/number.of.uniquely.mapped.bases"), what=numeric()))
}), use.names=TRUE)

combined.number.of.uniquely.mapped.bases.dt <- rbindlist(list(
    all.GSE36552.number.of.uniquely.mapped.bases.dt,
    all.GSE44183.number.of.uniquely.mapped.bases.dt
), use.names=TRUE)

fwrite(combined.number.of.uniquely.mapped.bases.dt, glue("{comparison.directory}/combined.number.of.uniquely.mapped.bases.count.dt.gz"))

{
    
    copy(combined.number.of.uniquely.mapped.bases.dt) -> .;
    .[, number.of.uniquely.mapped.bases.in.Gb:=(number.of.uniquely.mapped.bases/1e9)] -> ..to.plot.dt;
    .[, passed.threshold:=FALSE][number.of.uniquely.mapped.bases > 0.5*1e9, passed.threshold:=TRUE] -> .;
    ggplot(..to.plot.dt, aes(x=GSM, y=number.of.uniquely.mapped.bases.in.Gb, fill=passed.threshold)) -> .;
    . + geom_bar(stat="identity") -> .;
    . + geom_hline(yintercept=0.5, color="red", linetype="dashed") -> .;
    . + coord_flip() -> .;
    . + theme_pubr(base_size=10, legend="right") -> .;
    . + facet_grid(GSE~., scales="free", space="free") -> .;
    . + labs(x="GSM accession", y="# uniquely mapped bases in Gb\n(i.e., divided by 1E9)", fill="Passed the\nthreshold") -> .;
    ggsave.A4(filename=glue("{comparison.directory}/number.of.uniquely.mapped.bases.png"), plot=., width.r=0.45, height.r=0.8)

    
}


#### count of valid samples:
#### > nrow(combined.number.of.uniquely.mapped.bases.dt[number.of.uniquely.mapped.bases > 0.5 * 1000 * 1000 * 1000])
#### [1] 76

GSE36552.GSMs.passed.number.of.uniquely.mapped.bases.vector <- combined.number.of.uniquely.mapped.bases.dt[number.of.uniquely.mapped.bases > 0.5 * 1000 * 1000 * 1000][GSE == "GSE36552", GSM]
GSE44183.GSMs.passed.number.of.uniquely.mapped.bases.vector <- combined.number.of.uniquely.mapped.bases.dt[number.of.uniquely.mapped.bases > 0.5 * 1000 * 1000 * 1000][GSE == "GSE44183", GSM]

## 3. Get per-sample edits identified

all.valid.GSE36552.results.dt <- rbindlist(lapply(GSE36552.GSMs.passed.number.of.uniquely.mapped.bases.vector, function(temp.gsm){
    fread(glue("result/PS81_2__Qiu2016_filter_variants/single-100/201104-GSE36552-full124-100/{temp.gsm}/__merged__/trim-15bp-off-3prime-and-filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-2.8-1/GATK-2.8-1-default/151/common_all/dbSNP151.and.1000Genomeshg38.and.GoNLliftoverhg38/RNA_editing.sites.annotation.splicing.filter.txt"))[, `:=`(GSE="GSE36552", GSM=temp.gsm)]
}), use.names=TRUE)

all.valid.GSE44183.results.dt <- rbindlist(lapply(GSE44183.GSMs.passed.number.of.uniquely.mapped.bases.vector, function(temp.gsm){
    fread(glue("result/PS81_2__Qiu2016_filter_variants/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/{temp.gsm}/__merged__/filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-2.8-1/GATK-2.8-1-default/151/common_all/dbSNP151.and.1000Genomeshg38.and.GoNLliftoverhg38/RNA_editing.sites.annotation.splicing.filter.txt"))[, `:=`(GSE="GSE44183", GSM=temp.gsm)]
}), use.names=TRUE)


all.combined.valid.results.dt <- rbindlist(list(
    all.valid.GSE36552.results.dt,
    all.valid.GSE44183.results.dt
), use.names=TRUE)

fwrite(all.combined.valid.results.dt, glue("{comparison.directory}/all.combined.valid.results.dt.gz"))

#### median of edits across all samples
#### > median(all.combined.valid.results.dt[, .N, list(GSE, GSM)][, N])
#### [1] 3787.5

#### union of edits across all samples
#### > nrow(unique(all.combined.valid.results.dt[, list(Chromosome, Position)]))
#### [1] 271035


{

    copy(all.combined.valid.results.dt) -> .;
    .[, list(count=.N), list(GSE, GSM)] -> ..to.plot.dt;
    ggplot(..to.plot.dt, aes(x=GSM, y=count)) -> .;
    . + geom_bar(stat="identity") -> .;
    . + coord_flip() -> .;
    . + theme_pubr(base_size=10) -> .;
    . + facet_grid(GSE~., scales="free", space="free") -> .;
    . + labs(x="GSM accession", y="# edits identified") -> .;
    ggsave.A4(filename=glue("{comparison.directory}/count.of.edits.identified.raw.histogram.png"), plot=., width.r=0.45, height.r=0.85)

}


qiu2016.bed.liftOver.dt <- fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/qiu2016.bed.liftOver.directory/qiu2016.bed.lifOver.dt.gz")

{    
    merge(x=unique(all.combined.valid.results.dt[, list(Chromosome, Position, in.reimplementation=TRUE)]), y=qiu2016.bed.liftOver.dt[, list(Chr, pos, in.reported=TRUE)], by.x=c("Chromosome", "Position"), by.y=c("Chr", "pos"), all.x=TRUE, all.y=TRUE) -> .;
    .[is.na(in.reimplementation) == TRUE, in.reimplementation:=FALSE] -> .;
    .[is.na(in.reported) == TRUE, in.reported:=FALSE] -> .;
    .[, list(`detected in\nreimplementation`=in.reimplementation, `detected in\nQiu et al. reported`=in.reported)] -> ..to.plot.dt;
    png(glue("{comparison.directory}/comparison.between.reimplementation.and.Qiu2016.reported.png"), width=12, height=8, units="cm", res=600)
    print(plot(euler(..to.plot.dt), quantities=list(fontsize=7), legend=list(labels=c("detected in\nreimplementation", "detected in\nQiu et al. reported")), fills=c("lightblue", "#FF8F8F")))
    dev.off()   
}




#### overlap with the reported set
#### > merge(x=unique(all.combined.valid.results.dt[, list(Chromosome, Position, in.reimplementation=TRUE)]), y=qiu2016.bed.liftOver.dt[, list(Chr, pos, in.reported=TRUE)], by.x=c("Chromosome", "Position"), by.y=c("Chr", "pos"), all.x=TRUE, all.y=TRUE)[, .N, list(in.reimplementation, in.reported)]
####    in.reimplementation in.reported      N
#### 1:                TRUE          NA 267259
#### 2:                  NA        TRUE   5037
#### 3:                TRUE        TRUE   3776


## 4. Re-applying available filters and check again

#### 4.1. get the available filters first

all.GSE36552.MutDet.results.dt <- rbindlist(lapply(GSE36552.gsms.vector, function(temp.gsm){
    fread(glue("result/PS81_1__Qiu2016_call_variants/single-100/201104-GSE36552-full124-100/{temp.gsm}/__merged__/trim-15bp-off-3prime-and-filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-2.8-1/GATK-2.8-1-default/151/common_all/variation.sites.txt.gz"))[, `:=`(GSE="GSE36552", GSM=temp.gsm)]
}), use.names=TRUE)
setnames(all.GSE36552.MutDet.results.dt, c("Chr", "Pos", "Depth", "Ref", "Ref.Sup", "Ref.Mid", "Ref.End", "Ref.Plus", "Ref.Minus", "Ref.Mism", "Ref.NoMism", "Alt", "Alt.Sup", "Alt.Mid", "Alt.End", "Alt.Plus", "Alt.Minus", "Alt.Mism", "Alt.NoMism", "Index", "rposPvalue", "sbPvalue", "mmPvalue", "GSE", "GSM"))

all.GSE44183.MutDet.results.dt <- rbindlist(lapply(GSE44183.gsms.vector, function(temp.gsm){
    fread(glue("result/PS81_1__Qiu2016_call_variants/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/{temp.gsm}/__merged__/filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-2.8-1/GATK-2.8-1-default/151/common_all/variation.sites.txt.gz"))[, `:=`(GSE="GSE44183", GSM=temp.gsm)]
}), use.names=TRUE)
setnames(all.GSE44183.MutDet.results.dt, c("Chr", "Pos", "Depth", "Ref", "Ref.Sup", "Ref.Mid", "Ref.End", "Ref.Plus", "Ref.Minus", "Ref.Mism", "Ref.NoMism", "Alt", "Alt.Sup", "Alt.Mid", "Alt.End", "Alt.Plus", "Alt.Minus", "Alt.Mism", "Alt.NoMism", "Index", "rposPvalue", "sbPvalue", "mmPvalue", "GSE", "GSM"))

all.combined.MutDet.results.dt <- rbindlist(list(
    all.GSE36552.MutDet.results.dt,
    all.GSE44183.MutDet.results.dt
), use.names=TRUE)

all.combined.valid.results.with.Filter7.dt <- {
    copy(all.combined.valid.results.dt) -> .;
    ## metrics for Step 7
    .[, sample.occurrence:=length(unique(GSM)), list(Chromosome, Position, Reference, Alteration)] -> .;
    ## add Step7 filter
    .[, Step7.sample.occurrence.keep:=(sample.occurrence>=2)] -> .;
    ## return
    .
}

all.combined.valid.results.with.Filter12378.info.dt <- {
    copy(all.combined.MutDet.results.dt) -> .;
    ## metrics for Step 1 and Step 8
    .[, `:=`(Step.1.pvalue=mmPvalue, mapped.reads=Depth, variant.supporting.reads=Alt.Sup, mismatch.frequency=Alt.Sup/Depth)] -> .;
    ## metrics for Step 2
    .[, `:=`(Step.2.pvalue=sbPvalue, variant.strand.frequency=Alt.Plus/Alt.Sup, reference.strand.frequency=Ref.Plus/Ref.Sup)] -> .;
    .[, `:=`(variant.strand.preference=abs(variant.strand.frequency - 0.5), reference.strand.preference=abs(reference.strand.frequency - 0.5))] -> .;
    ## metrics for Step 3
    .[, `:=`(Step.3.pvalue=rposPvalue, read.end.frequency=Alt.End/Alt.Sup, read.middle.frequency=Alt.Mid/Alt.Sup)] -> .;
    ## determine the filters
    .[,`:=`(
        Step1.pvalue.keep=mmPvalue<0.01,
        Step1.mapped.reads.keep=mapped.reads>=4,
        Step1.variant.supporting.reads.keep=variant.supporting.reads>=3,
        Step1.mismatch.frequency.keep=mismatch.frequency>=0.1,
        Step2.pvalue.keep=((sbPvalue<0.005) & (variant.strand.preference > reference.strand.preference)) == FALSE,
        Step2.variant.strand.frequency.not.too.high.keep=((variant.strand.frequency > 0.9) == FALSE),
        Step2.variant.strand.frequency.not.too.low.keep=((variant.strand.frequency < 0.1) == FALSE),
        Step3.pvalue.keep=((rposPvalue<0.05) & (read.end.frequency > read.middle.frequency)) == FALSE,
        Step3.read.end.frequency.keep=(read.end.frequency > 0.9) == FALSE,
        Step8.mismatch.frequency.not.too.high.keep=(mismatch.frequency>0.95)==FALSE
    )] -> .;
    ## LEFT-JOIN all.combined.valid.results.dt WITH .
    merge(x=all.combined.valid.results.with.Filter7.dt, y=.,
          by.x=c("GSE", "GSM", "Chromosome", "Position"),
          by.y=c("GSE", "GSM", "Chr", "Pos"),
          all.x=TRUE, all.y=FALSE) -> .;
    .
}

#### 4.2. apply the filters

all.combined.valid.results.with.Filter12378.info.refiltered.dt <- all.combined.valid.results.with.Filter12378.info.dt[(
    Step1.pvalue.keep & Step1.mapped.reads.keep & Step1.variant.supporting.reads.keep & Step1.mismatch.frequency.keep &
    Step2.pvalue.keep & Step2.variant.strand.frequency.not.too.high.keep & Step2.variant.strand.frequency.not.too.low.keep &
    Step3.pvalue.keep & Step3.read.end.frequency.keep &
    Step7.sample.occurrence.keep &
    Step8.mismatch.frequency.not.too.high.keep
)]

fwrite(all.combined.valid.results.with.Filter12378.info.refiltered.dt, glue("{comparison.directory}/all.combined.valid.results.with.Filter12378.info.refiltered.dt.gz"))

{

    copy(all.combined.valid.results.with.Filter12378.info.refiltered.dt) -> .;
    .[, list(count=.N), list(GSE, GSM)] -> ..to.plot.dt;
    ggplot(..to.plot.dt, aes(x=GSM, y=count)) -> .;
    . + geom_bar(stat="identity") -> .;
    . + coord_flip() -> .;
    . + theme_pubr(base_size=10) -> .;
    . + facet_grid(GSE~., scales="free", space="free") -> .;
    . + labs(x="GSM accession", y="# edits identified after re-filtering") -> .;
    ggsave.A4(filename=glue("{comparison.directory}/count.of.edits.identified.refiltered.histogram.png"), plot=., width.r=0.45, height.r=0.85)

}


{    
    merge(x=unique(all.combined.valid.results.with.Filter12378.info.refiltered.dt[, list(Chromosome, Position, in.reimplementation=TRUE)]), y=qiu2016.bed.liftOver.dt[, list(Chr, pos, in.reported=TRUE)], by.x=c("Chromosome", "Position"), by.y=c("Chr", "pos"), all.x=TRUE, all.y=TRUE) -> .;
    .[is.na(in.reimplementation) == TRUE, in.reimplementation:=FALSE] -> .;
    .[is.na(in.reported) == TRUE, in.reported:=FALSE] -> .;
    .[, list(`detected in\nreimplementation (re-filtered)`=in.reimplementation, `detected in\nQiu et al. reported`=in.reported)] -> ..to.plot.dt;
    png(glue("{comparison.directory}/comparison.between.reimplementation.refiltered.and.Qiu2016.reported.png"), width=12, height=8, units="cm", res=600)
    print(plot(euler(..to.plot.dt), quantities=list(fontsize=7), legend=list(labels=c("detected in\nreimplementation\n(re-filtered)", "detected in\nQiu et al. reported")), fills=c("lightblue", "#FF8F8F")))
    dev.off()   
}


#### median of edits across all samples
#### > median(all.combined.valid.results.with.Filter12378.info.refiltered.dt[, .N, list(GSE, GSM)][, N])
#### [1] 49

#### union of edits across all samples
#### > nrow(unique(all.combined.valid.results.with.Filter12378.info.refiltered.dt[, list(Chromosome, Position)]))
#### [1] 2700

#### overlap with the reported set
#### > merge(x=unique(all.combined.valid.results.with.Filter12378.info.refiltered.dt[, list(Chromosome, Position, in.reimplementation=TRUE)]), y=qiu2016.bed.liftOver.dt[, list(Chr, pos, in.reported=TRUE)], by.x=c("Chromosome", "Position"), by.y=c("Chr", "pos"), all.x=TRUE, all.y=TRUE)[, .N, list(in.reimplementation, in.reported)]
####    in.reimplementation in.reported    N
#### 1:                TRUE          NA 1804
#### 2:                  NA        TRUE 7917
#### 3:                TRUE        TRUE  896



## 5. Determine the final set of Qiu et al.'s edits

Qiu2016.edits.dt <- all.combined.valid.results.with.Filter12378.info.dt[(
    Step1.mapped.reads.keep & Step1.variant.supporting.reads.keep & Step1.mismatch.frequency.keep &
    Step2.pvalue.keep & Step2.variant.strand.frequency.not.too.high.keep & Step2.variant.strand.frequency.not.too.low.keep &
    Step3.pvalue.keep & Step3.read.end.frequency.keep &
    Step7.sample.occurrence.keep &
    Step8.mismatch.frequency.not.too.high.keep
)]

fwrite(Qiu2016.edits.dt, glue("{comparison.directory}/Qiu2016.edits.dt.gz"))

{

    copy(Qiu2016.edits.dt) -> .;
    .[, list(count=.N), list(GSE, GSM)] -> ..to.plot.dt;
    ggplot(..to.plot.dt, aes(x=GSM, y=count)) -> .;
    . + geom_bar(stat="identity") -> .;
    . + coord_flip() -> .;
    . + theme_pubr(base_size=10) -> .;
    . + facet_grid(GSE~., scales="free", space="free") -> .;
    . + labs(x="GSM accession", y="# edits identified (final)") -> .;
    ggsave.A4(filename=glue("{comparison.directory}/count.of.edits.identified.finalset.histogram.png"), plot=., width.r=0.45, height.r=0.85)

}


{    
    merge(x=unique(Qiu2016.edits.dt[, list(Chromosome, Position, in.reimplementation=TRUE)]), y=qiu2016.bed.liftOver.dt[, list(Chr, pos, in.reported=TRUE)], by.x=c("Chromosome", "Position"), by.y=c("Chr", "pos"), all.x=TRUE, all.y=TRUE) -> .;
    .[is.na(in.reimplementation) == TRUE, in.reimplementation:=FALSE] -> .;
    .[is.na(in.reported) == TRUE, in.reported:=FALSE] -> .;
    .[, list(`detected in\nreimplementation (final)`=in.reimplementation, `detected in\nQiu et al. reported`=in.reported)] -> ..to.plot.dt;
    png(glue("{comparison.directory}/comparison.between.reimplementation.finalset.and.Qiu2016.reported.png"), width=12, height=8, units="cm", res=600)
    print(plot(euler(..to.plot.dt), quantities=list(fontsize=7), legend=list(labels=c("detected in\nreimplementation\n(final)", "detected in\nQiu et al. reported")), fills=c("lightblue", "#FF8F8F")))
    dev.off()   
}

#### median of edits across all samples
#### > median(Qiu2016.edits.dt[, .N, list(GSE, GSM)][, N])
#### [1] 449.5

#### union of edits across all samples
#### > nrow(unique(Qiu2016.edits.dt[, list(Chromosome, Position)]))
#### [1] 13500

#### overlap with the reported set
#### > merge(x=unique(Qiu2016.edits.dt[, list(Chromosome, Position, in.reimplementation=TRUE)]), y=qiu2016.bed.liftOver.dt[, list(Chr, pos, in.reported=TRUE)], by.x=c("Chromosome", "Position"), by.y=c("Chr", "pos"), all.x=TRUE, all.y=TRUE)[, .N, list(in.reimplementation, in.reported)]
####    in.reimplementation in.reported     N
#### 1:                TRUE          NA 10165
#### 2:                  NA        TRUE  5478
#### 3:                TRUE        TRUE  3335


## 6. Compare with ours and find out FN and FP sets (assuming that Qiu2016 is the real dataset)

our.edits.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")
our.edits.GSE36552.and.GSE44183.dt <- our.edits.dt[SAMPLE %in% c(GSE36552.GSMs.passed.number.of.uniquely.mapped.bases.vector, GSE44183.GSMs.passed.number.of.uniquely.mapped.bases.vector)]

Qiu2016.edits.and.ours.dt <- {
    copy(Qiu2016.edits.dt) -> .;
    .[, detected.in.Qiu2016:=TRUE] -> .;
    merge(x=., y=data.table(our.edits.GSE36552.and.GSE44183.dt, detected.in.ours=TRUE),
          by.x=c("GSE", "GSM", "Chromosome", "Position"),
          by.y=c("gse", "SAMPLE", "CHROM", "POS"),
          all.x=TRUE, all.y=TRUE) -> .;
    .[is.na(detected.in.Qiu2016) == TRUE, detected.in.Qiu2016:=FALSE] -> .;
    .[is.na(detected.in.ours) == TRUE, detected.in.ours:=FALSE] -> .;
    .    
}

fwrite(Qiu2016.edits.and.ours.dt, glue("{comparison.directory}/Qiu2016.edits.and.ours.dt.gz"))

{

    copy(Qiu2016.edits.and.ours.dt) -> .;
    .[, list(count=.N), list(GSE, GSM, detected.in.ours, detected.in.Qiu2016)] -> .;
    .[detected.in.ours == TRUE & detected.in.Qiu2016 == TRUE, group:="TP: detected\nin both"] -> .;
    .[detected.in.ours == FALSE & detected.in.Qiu2016 == TRUE, group:="FN: detected\nin Qiu et al. only"] -> .;
    .[detected.in.ours == TRUE & detected.in.Qiu2016 == FALSE, group:="FP: detected\nin ours only"] -> .;
    . -> ..to.plot.dt

    ggplot(..to.plot.dt, aes(x=GSM, y=count, fill=group)) -> .;
    . + geom_bar(stat="identity", position="fill") -> .;
    . + coord_flip() -> .;
    . + theme_pubr(base_size=10) -> .;
    . + facet_grid(GSE~., scales="free", space="free") -> .;
    . + labs(x="GSM accession", y="Ratio") -> .;
    . + scale_fill_manual(values=c("#8F8FFF", "#FF8F8F", "#8FFF8F")) -> .;
    ggsave.A4(filename=glue("{comparison.directory}/overlap.of.edits.identified.between.finalset.and.ours.histogram.png"), plot=., width.r=0.75, height.r=0.85)

}


#### summarize |intersection = detected in both|/|union = detected in at least one set| for each sample
#### > Qiu2016.edits.and.ours.dt[, list(count=.N), list(detected.in.ours, detected.in.Qiu2016, GSE, GSM)][, total.count.per.GSM:=sum(count), list(GSM)][, percentage:=count/total.count.per.GSM*100][detected.in.Qiu2016 == TRUE & detected.in.ours == TRUE, data.table(t(quantile(percentage)))]
####           0%     25%     50%      75%     100%
#### 1: 0.6802721 4.53691 5.85434 8.408017 10.29439



## 7. Check false negatives:

ours.false.negatives.dt <- Qiu2016.edits.and.ours.dt[detected.in.ours == FALSE]
ours.false.negatives.dt[, Chromosome.and.Position:=paste(sep="", Chromosome, "_", Position)]

#### check total and median
#### > ours.false.negatives.dt[, length(unique(Chromosome.and.Position))]
#### [1] 9798
#### > ours.false.negatives.dt[, .N, GSM][, median(N)]
#### [1] 251



## Check FN-1: overlapping with cohort
all.chromosomes.vector <- glue("chr{c(1:22, 'X', 'Y')}")
for(chr in all.chromosomes.vector){
    fwrite(unique(ours.false.negatives.dt[Chromosome == chr, list(Chromosome, start=Position-1, end=Position-1+1)]), glue("{comparison.directory}/ours.false.negatives.unique.location.{chr}.bed"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

registerDoMC(cores=10)
foreach(temp.chr=all.chromosomes.vector) %dopar% {
    temp.outers.collection <- c("1000Genomes.phased.genotype", "gnomAD_v2.1.1_exomes", "gnomAD_v2.1.1_genomes", "gnomAD_v3.0_genomes", "NCBI.ALFA.2020.03.04", "UWashington.EVS")
    if (temp.chr == "chrY") {
        temp.outers.collection <- c("gnomAD_v2.1.1_exomes", "gnomAD_v3.0_genomes", "NCBI.ALFA.2020.03.04", "UWashington.EVS")
    }
    foreach(temp.outer=temp.outers.collection) %do% {
        temp.comparison.directory <- glue("{comparison.directory}/outer.vcf.comparison/{temp.outer}/{temp.chr}/")
        dir.create(temp.comparison.directory, recursive=TRUE)
        system(glue("~/tool/bcftools-1.15/bcftools view --regions-file {comparison.directory}/ours.false.negatives.unique.location.{temp.chr}.bed -o {temp.comparison.directory}/temp.out.vcf ./external/outer_vcf/{temp.outer}/{temp.chr}/outer.VCF"))
    }
} -> not.used.variable;

rbindlist(foreach(temp.chr=all.chromosomes.vector) %do% {
    temp.outers.collection <- c("1000Genomes.phased.genotype", "gnomAD_v2.1.1_exomes", "gnomAD_v2.1.1_genomes", "gnomAD_v3.0_genomes", "NCBI.ALFA.2020.03.04", "UWashington.EVS")
    if (temp.chr == "chrY") {
        temp.outers.collection <- c("gnomAD_v2.1.1_exomes", "gnomAD_v3.0_genomes", "NCBI.ALFA.2020.03.04", "UWashington.EVS")
    }
    rbindlist(foreach(temp.outer=temp.outers.collection) %do% {
        temp.comparison.directory <- glue("{comparison.directory}/outer.vcf.comparison/{temp.outer}/{temp.chr}/")
        dir.create(temp.comparison.directory, recursive=TRUE)
        fread(cmd=glue("grep -v '^##' {temp.comparison.directory}/temp.out.vcf | cut -f 1-3"))[, outer:=temp.outer]
    }, use.names=TRUE)
}, use.names=TRUE) -> outer.vcf.overlap.with.ours.false.negatives.dt;
outer.vcf.overlap.with.ours.false.negatives.dt[, Chr.and.Pos:=paste(sep="", `#CHROM`, "_", POS)]


ours.false.negatives.dt[, overlaps.with.outer.vcf:=FALSE][Chromosome.and.Position %in% outer.vcf.overlap.with.ours.false.negatives.dt[, Chr.and.Pos], overlaps.with.outer.vcf:=TRUE]

## Check FN-2: multiple alternative alleles


multiple.alternative.allele.sites.in.ours.dt <- fread(cmd="zcat result/S51_2__filter_against_population_variants/210215-sixth-dataset/merged.variant.only.disjoint.with.population.variants.vcf.gz | grep -v '^##' | awk '{if ($3 ~ /;/) {print $1\"\t\"$2\"\t\"$3}}'")
setnames(multiple.alternative.allele.sites.in.ours.dt, c("Chr", "Pos", "ID"))
multiple.alternative.allele.sites.in.ours.dt[, Chr.and.Pos:=paste(sep="", Chr, "_", Pos)]

ours.false.negatives.dt[, overlaps.with.multiple.alternative.allele.sites:=FALSE][Chromosome.and.Position %in% multiple.alternative.allele.sites.in.ours.dt[, Chr.and.Pos], overlaps.with.multiple.alternative.allele.sites:=TRUE]




## Check FN-3: incorrect strands

combined.merged.variant.only.snpEff.event.summary.dt <- fread("result/S18_1__combine_annotations/201218-fifth-dataset/__merged__/base-quality-no-smaller-than-25/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/snpEff/basic/10000000/combined.merged.variant.only.snpEff.event.summary.dt.txt.gz")
combined.merged.variant.only.snpEff.event.summary.dt[, Chr.and.Pos:=paste(sep="", CHROM, "_", POS)]
combined.merged.variant.only.snpEff.event.summary.A.to.G.dt <- combined.merged.variant.only.snpEff.event.summary.dt[event.summary %in% c("A>G", "A>G;T>C")]

ours.false.negatives.dt[, event.does.not.contain.A.to.G:=TRUE][Chromosome.and.Position %in% combined.merged.variant.only.snpEff.event.summary.A.to.G.dt[, Chr.and.Pos], event.does.not.contain.A.to.G:=FALSE]

## Check FN-4 : sample occurrence

merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt <- fread("result/S51_4__filter_for_variants_with_enough_sample_support/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt.txt.gz")
sample.occurrence.per.edits.for.GSE36552.and.GSE44183.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt[gse %in% c("GSE36552", "GSE44183")][, Chr.and.Pos:=sub(pattern="^(chr[0-9XY]+_[0-9]+)_.*", replacement="\\1", x=ID)]

ours.false.negatives.dt[, invalid.sample.occurrence:=TRUE][Chromosome.and.Position %in% sample.occurrence.per.edits.for.GSE36552.and.GSE44183.dt[(SUBSET =="Alu") | (SUBSET != "Alu" & site.occurrence >=2), Chr.and.Pos], invalid.sample.occurrence:=FALSE]

## write the results

fwrite(ours.false.negatives.dt, glue("{comparison.directory}/ours.false.negatives.dt.gz"))

{

    copy(ours.false.negatives.dt) -> .;
    ## updated order: FN-2, FN-1, FN-3, FN-4
    .[, reason:="5. Others"] -> .;
    .[overlaps.with.multiple.alternative.allele.sites == TRUE, reason:="1. Multiple alleles detected"] -> .;
    .[overlaps.with.multiple.alternative.allele.sites == FALSE & overlaps.with.outer.vcf == TRUE, reason:="2. Overlapping with variants\nfrom worldwide cohort genotyping"] -> .;
    .[overlaps.with.outer.vcf == FALSE & overlaps.with.multiple.alternative.allele.sites == FALSE & invalid.sample.occurrence == TRUE, reason:="3. Invalid sample\noccurrence"] -> .;
    .[overlaps.with.outer.vcf == FALSE & overlaps.with.multiple.alternative.allele.sites == FALSE & invalid.sample.occurrence == FALSE & event.does.not.contain.A.to.G == TRUE, reason:="4. No A-to-G variant\non known transcripts"] -> .;
    .[, list(count=.N), list(GSE, GSM, reason)][, total.count:=sum(count), list(GSM)][, percentage:=count/total.count*100] -> .;
    . -> ..to.plot.dt
    ##
    ggplot(..to.plot.dt, aes(x=GSM, y=count, fill=reason)) -> .;
    . + geom_bar(stat="identity", position="fill") -> .;
    . + coord_flip() -> .;
    . + theme_pubr(base_size=10) -> .;
    . + geom_hline(yintercept=0.25, color="black", linetype="dashed") -> .;
    . + guides(fill=guide_legend(nrow=2)) -> .;
    . + facet_grid(GSE~., scales="free", space="free") -> .;
    . + labs(x="GSM accession", y="Ratio", fill="Reasons of FN") -> .;
    . + scale_fill_manual(values=c("#0000FF", "#5050FF", "#A0A0FF", "#F0F0FF", "grey50")) -> .;
    ggsave.A4(filename=glue("{comparison.directory}/false.negative.explanation.png"), plot=., width.r=0.75, height.r=0.85)

}



## ## summarization
## ours.false.negatives.summary.dt <- {
##     copy(ours.false.negatives.dt) -> .;
##     summarize the per-step explanation
##     .[, reason:="5. Others"] -> .;
##     .[overlaps.with.multiple.alternative.allele.sites == TRUE, reason:="1. Multiple alleles detected"] -> .;
##     .[overlaps.with.multiple.alternative.allele.sites == FALSE & overlaps.with.outer.vcf == TRUE, reason:="2. Overlapping with variants\nfrom worldwide cohort genotyping"] -> .;
##     .[overlaps.with.outer.vcf == FALSE & overlaps.with.multiple.alternative.allele.sites == FALSE & invalid.sample.occurrence == TRUE, reason:="3. Invalid sample\noccurrence"] -> .;
##     .[overlaps.with.outer.vcf == FALSE & overlaps.with.multiple.alternative.allele.sites == FALSE & invalid.sample.occurrence == FALSE & event.does.not.contain.A.to.G == TRUE, reason:="4. No A-to-G variant\non known transcripts"] -> .;
##     .[, list(count=.N), list(GSE, GSM, reason)][, total.count:=sum(count), list(GSM)][, percentage:=count/total.count*100] -> .;
##     .
##  }

#### check reason distribution
#### > ours.false.negatives.summary.dt[, data.table(t(quantile(percentage))), reason][order(reason)]
##                                                            reason        0%       25%      50%       75%      100%
## 1:                                   1. Multiple alleles detected  1.562500  3.298611  3.90625  4.794521  7.751938
## 2: 2. Overlapping with variants\nfrom worldwide cohort genotyping 41.787440 49.255023 54.57589 65.632393 76.635514
## 3:                                  3. Invalid sample\noccurrence  2.419355  7.928604 10.13294 11.923155 29.411765
## 4:                     4. No A-to-G variant\non known transcripts  1.851852 10.756501 20.28189 23.616863 31.707317
## 5:                                                      5. Others  3.174603  8.935830 12.10191 16.187439 24.050633




## 8: Check false positives

ours.false.positives.dt <- Qiu2016.edits.and.ours.dt[detected.in.Qiu2016 == FALSE]
ours.false.positives.dt[, Chromosome.and.Position:=paste(sep="", Chromosome, "_", Position)]

#### check total and median
#### > ours.false.positives.dt[, length(unique(Chromosome.and.Position))]
#### [1] 108877
#### > ours.false.positives.dt[, .N, GSM][, median(N)]
#### [1] 3058


## Check FP-1: reads mappable by bwa but unmappable by tophat-2


registerDoMC(cores=6)
## generating pileup files
foreach(temp.row.dt=iter(GSE.and.GSM.dt[GSM %in% c(GSE36552.GSMs.passed.number.of.uniquely.mapped.bases.vector, GSE44183.GSMs.passed.number.of.uniquely.mapped.bases.vector)], by="row")) %dopar% {
    temp.GSE <- temp.row.dt[1, GSE]
    temp.GSM <- temp.row.dt[1, GSM]
    cat(glue("{date()} Processing GSE {temp.GSE} and GSM {temp.GSM} ..."), "\n")
    temp.output.directory <- glue("{comparison.directory}/samtools.mpileup.comparison/{temp.GSE}/{temp.GSM}")
    dir.create(temp.output.directory, recursive=TRUE)
    temp.FP.sites.bed.filename <- glue("{temp.output.directory}/FP.sites.bed")
    fwrite(ours.false.positives.dt[GSE==temp.GSE & GSM==temp.GSM, list(Chromosome, start=Position-1, end=Position-1+1)], temp.FP.sites.bed.filename, row.names=FALSE, col.names=FALSE, sep="\t")
    ## set up names
    temp.ours.alignment.bam <- NA
    temp.ours.GATK.recal.bam <- NA
    temp.Qiu2016.alignment.bam <- NA
    temp.Qiu2016.bwa.realignment.readS.bam <- NA
    temp.Qiu2016.bwa.realignment.readP.bam <- NA
    if (temp.GSE == "GSE36552") {
        temp.ours.alignment.bam <- glue("result/S15_1__get_sample_RNA_editing_sites_v3/single-100/201104-GSE36552-full124-100/{temp.GSM}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/95/bwa-aln-samsepe/none/alignment.merged.bam")
        temp.ours.GATK.recal.bam <- glue("result/S15_1__get_sample_RNA_editing_sites_v3/single-100/201104-GSE36552-full124-100/{temp.GSM}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/95/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam")
        temp.Qiu2016.alignment.bam <- glue("result/PS15_1__get_sample_RNA_editing_sites_v3/single-100/201104-GSE36552-full124-100/{temp.GSM}/__merged__/trim-15bp-off-3prime-and-filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/alignment.merged.bam")
        temp.Qiu2016.bwa.realignment.readS.bam <- glue("result/PS81_2__Qiu2016_filter_variants/single-100/201104-GSE36552-full124-100/{temp.GSM}/__merged__/trim-15bp-off-3prime-and-filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-2.8-1/GATK-2.8-1-default/151/common_all/dbSNP151.and.1000Genomeshg38.and.GoNLliftoverhg38/bwa/mutation.readS.bam")
    } else if (temp.GSE == "GSE44183") {
        temp.ours.alignment.bam <- glue("result/S15_1__get_sample_RNA_editing_sites_v3/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/{temp.GSM}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/85/bwa-aln-samsepe/none/alignment.merged.bam")
        temp.ours.GATK.recal.bam <- glue("result/S15_1__get_sample_RNA_editing_sites_v3/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/{temp.GSM}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/85/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam")
        temp.Qiu2016.alignment.bam <- glue("result/PS15_1__get_sample_RNA_editing_sites_v3/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/{temp.GSM}/__merged__/filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/alignment.merged.bam")
        temp.Qiu2016.bwa.realignment.readS.bam <- glue("result/PS81_2__Qiu2016_filter_variants/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/{temp.GSM}/__merged__/filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-2.8-1/GATK-2.8-1-default/151/common_all/dbSNP151.and.1000Genomeshg38.and.GoNLliftoverhg38/bwa/mutation.readS.bam")
        temp.Qiu2016.bwa.realignment.readP.bam <- glue("result/PS81_2__Qiu2016_filter_variants/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/{temp.GSM}/__merged__/filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-2.8-1/GATK-2.8-1-default/151/common_all/dbSNP151.and.1000Genomeshg38.and.GoNLliftoverhg38/bwa/mutation.readP.bam")
    }
    ## check raw bam
    system(glue("samtools mpileup --output-QNAME -l {temp.FP.sites.bed.filename} {temp.ours.alignment.bam} > {temp.output.directory}/FP.sites.ours.mpileup"))
    system(glue("samtools mpileup --output-QNAME -l {temp.FP.sites.bed.filename} {temp.Qiu2016.alignment.bam} > {temp.output.directory}/FP.sites.Qiu2016.mpileup"))
    ## check our GATK recal bam
    system(glue("samtools mpileup --output-QNAME -l {temp.FP.sites.bed.filename} {temp.ours.GATK.recal.bam} > {temp.output.directory}/FP.sites.ours.GATK.recal.mpileup")) 
    ## check bwa realignment readS
    system(glue("samtools sort -o {temp.output.directory}/temp.Qiu2016.bwa.realignment.readS.sorted.bam {temp.Qiu2016.bwa.realignment.readS.bam}"))
    system(glue("samtools mpileup --output-QNAME -l {temp.FP.sites.bed.filename} {temp.output.directory}/temp.Qiu2016.bwa.realignment.readS.sorted.bam > {temp.output.directory}/FP.sites.Qiu2016.bwa.realignment.readS.mpileup"))
    ## (GSE44183 only) check bwa realignment readP
    if (temp.GSE == "GSE44183") {
        system(glue("samtools sort -o {temp.output.directory}/temp.Qiu2016.bwa.realignment.readP.sorted.bam {temp.Qiu2016.bwa.realignment.readP.bam}"))
        system(glue("samtools mpileup --output-QNAME -l {temp.FP.sites.bed.filename} {temp.output.directory}/temp.Qiu2016.bwa.realignment.readP.sorted.bam > {temp.output.directory}/FP.sites.Qiu2016.bwa.realignment.readP.mpileup"))
    }
} -> not.used.variable



combined.mpileup.dt <- rbindlist(use.names=TRUE, foreach(temp.row.dt=iter(GSE.and.GSM.dt[GSM %in% c(GSE36552.GSMs.passed.number.of.uniquely.mapped.bases.vector, GSE44183.GSMs.passed.number.of.uniquely.mapped.bases.vector)], by="row")) %do% {
    temp.GSE <- temp.row.dt[1, GSE]
    temp.GSM <- temp.row.dt[1, GSM]
    cat(glue("{date()} Reading GSE {temp.GSE} and GSM {temp.GSM} ..."), "\n")
    temp.output.directory <- glue("{comparison.directory}/samtools.mpileup.comparison/{temp.GSE}/{temp.GSM}")
    rbindlist(list(
        fread(glue("{temp.output.directory}/FP.sites.ours.mpileup"), header=FALSE)[, group:="ours"],
        fread(glue("{temp.output.directory}/FP.sites.ours.GATK.recal.mpileup", header=FALSE))[, group:="ours.GATK.recal"],
        fread(glue("{temp.output.directory}/FP.sites.Qiu2016.mpileup"), header=FALSE)[, group:="Qiu2016"],
        fread(glue("{temp.output.directory}/FP.sites.Qiu2016.bwa.realignment.readS.mpileup"), header=FALSE)[, group:="Qiu2016.bwa.realignment.readS"],
        {if (temp.GSE == "GSE44183") {fread(glue("{temp.output.directory}/FP.sites.Qiu2016.bwa.realignment.readP.mpileup"), header=FALSE)[, group:="Qiu2016.bwa.realignment.readP"]} else {NULL} }
    ), use.names=TRUE) -> temp.mpileup.dt
    setnames(temp.mpileup.dt, c("Chromosome", "Position", "Ref", "Depth", "Bases", "Base.Qualities", "QNAMEs", "group"))
    temp.mpileup.dt[, `:=`(GSE=temp.GSE, GSM=temp.GSM)]
    temp.mpileup.dt
})
combined.mpileup.dt[, read.count:=str_count(string=QNAMEs, pattern=",") + 1]
combined.mpileup.dcast.dt <- dcast(combined.mpileup.dt, GSE + GSM + Chromosome + Position ~ group, value.var="read.count", fill=0)
combined.mpileup.dcast.dt[, `:=`(Qiu2016.not.detected.in.first.alignment=(Qiu2016==0), Qiu2016.not.detected.in.bwa.realignment=(Qiu2016.bwa.realignment.readP == 0 & Qiu2016.bwa.realignment.readS==0))]

combined.Qiu2016.bwa.filtered.dt <- rbindlist(use.names=TRUE, foreach(temp.row.dt=iter(GSE.and.GSM.dt[GSM %in% c(GSE36552.GSMs.passed.number.of.uniquely.mapped.bases.vector, GSE44183.GSMs.passed.number.of.uniquely.mapped.bases.vector)], by="row")) %do% {
    temp.GSE <- temp.row.dt[1, GSE]
    temp.GSM <- temp.row.dt[1, GSM]
    cat(glue("{date()} Processing GSE {temp.GSE} and GSM {temp.GSM} ..."), "\n")
    if (temp.GSE=="GSE36552"){
        fread(glue("result/PS81_2__Qiu2016_filter_variants/single-100/201104-GSE36552-full124-100/{temp.GSM}/__merged__/trim-15bp-off-3prime-and-filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-2.8-1/GATK-2.8-1-default/151/common_all/dbSNP151.and.1000Genomeshg38.and.GoNLliftoverhg38/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.txt"), header=FALSE, select=c(1,2), col.names=c("Chr", "Pos"))[, `:=`(GSE=temp.GSE, GSM=temp.GSM)]
    } else if (temp.GSE=="GSE44183"){
        fread(glue("result/PS81_2__Qiu2016_filter_variants/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/{temp.GSM}/__merged__/filter-by-soapnuke/hg38.fa/32/tophat2-index/tophat2-index-default/tophat2/tophat2-default/GATK-2.8-1/GATK-2.8-1-default/151/common_all/dbSNP151.and.1000Genomeshg38.and.GoNLliftoverhg38/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.txt"), header=FALSE, select=c(1,2), col.names=c("Chr", "Pos"))[, `:=`(GSE=temp.GSE, GSM=temp.GSM)]
    } else{
        stop(glue("Unsupported GSE {temp.GSE}"))
    }
})
combined.mpileup.dcast.dt[, discarded.in.bwa:=TRUE][paste(sep="_", GSE, GSM, Chromosome, Position) %in% combined.Qiu2016.bwa.filtered.dt[, paste(sep="_", GSE, GSM, Chr, Pos)], discarded.in.bwa:=FALSE]

fwrite(combined.mpileup.dcast.dt, glue("{comparison.directory}/combined.mpileup.dcast.dt.gz"))


{

    copy(combined.mpileup.dcast.dt) -> .;
    .[, reason:="4. Others"] ->.;
    .[Qiu2016.not.detected.in.first.alignment == TRUE, reason:="1. Not detected in first-round alignment"] -> .;
    .[Qiu2016.not.detected.in.first.alignment == FALSE & Qiu2016.not.detected.in.bwa.realignment == TRUE, reason:="2. Not detected in BWA realignment"] -> .;
    .[Qiu2016.not.detected.in.first.alignment == FALSE & Qiu2016.not.detected.in.bwa.realignment == FALSE & discarded.in.bwa==TRUE, reason:="3. Detected in BWA realignment\nbut discarded anyway"] -> .;
    .[, list(count=.N), list(GSE, GSM, reason)][, total.count:=sum(count), list(GSM)][, percentage:=count/total.count*100] -> .;
    . -> ..to.plot.dt

    ggplot(..to.plot.dt, aes(x=GSM, y=count, fill=reason)) -> .;
    . + geom_bar(stat="identity", position="fill") -> .;
    . + coord_flip() -> .;
    . + theme_pubr(base_size=10) -> .;
    . + geom_hline(yintercept=0.34, color="black", linetype="dashed") -> .;
    . + guides(fill=guide_legend(nrow=2)) -> .;
    . + facet_grid(GSE~., scales="free", space="free") -> .;
    . + labs(x="GSM accession", y="Ratio", fill="Reasons of FP") -> .;
    . + scale_fill_manual(values=c("#FF0000", "#FF5A5A", "#FFB4B4", "grey50")) -> .;
    ggsave.A4(filename=glue("{comparison.directory}/false.positive.explanation.png"), plot=., width.r=0.75, height.r=0.85)

}



## combined.mpileup.dcast.summary.dt <- {
##     copy(combined.mpileup.dcast.dt) -> .;
##     .[, reason:="4.others"] ->.;
##     .[Qiu2016.not.detected.in.first.alignment == TRUE, reason:="1.Qiu2016.not.detected.in.first.alignment"] -> .;
##     .[Qiu2016.not.detected.in.first.alignment == FALSE & Qiu2016.not.detected.in.bwa.realignment == TRUE, reason:="2.Qiu2016.not.detected.in.bwa.realignment"] -> .;
##     .[Qiu2016.not.detected.in.first.alignment == FALSE & Qiu2016.not.detected.in.bwa.realignment == FALSE & discarded.in.bwa==TRUE, reason:="3.Qiu2016.detected.but.discarded.in.bwa.realignment"] -> .;
##     .
## }

#### summarize reasons

#### > combined.mpileup.dcast.summary.dt[, list(count=.N), list(GSE, GSM, reason)][, total.count:=sum(count), list(GSM)][, percentage:=count/total.count * 100][, data.table(t(quantile(percentage))), list(reason)][order(reason)]
####                                                 reason         0%       25%       50%       75%     100%
#### 1:           1.Qiu2016.not.detected.in.first.alignment  0.1494449  0.593958  2.806019  6.228692 53.15568
#### 2:           2.Qiu2016.not.detected.in.bwa.realignment 16.5497896 29.270797 34.964236 42.007955 55.17241
#### 3: 3.Qiu2016.detected.but.discarded.in.bwa.realignment 21.0640608 29.018609 40.482429 49.164659 58.54813
#### 4:                                            4.others  1.3537906  5.419275 17.218312 33.431319 45.27687

