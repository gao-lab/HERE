library("data.table")
library("glue")

"GSM2706237" -> sample.name
"result/S15_1__get_sample_RNA_editing_sites_v3/paired-100-100/200902-GSE101571-full-100-100/GSM2706237/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/95/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all" -> sample.specific.S15.path
"result/S15_1__get_sample_RNA_editing_sites_v3/paired-100-100/200902-GSE101571-full-100-100/GSM2706237/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/95/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none" -> sample.specific.S15.deeper.path

{

    ## raw.pct
    glue("bcftools view {sample.specific.S15.path}/alignment.bcf |grep -v '^##' ") -> .;
    fread(cmd=., header=TRUE) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> raw.pct
    ## print(raw.pct)

    ## after step 1
    glue("{sample.specific.S15.deeper.path}/alignment.con.vcf") -> .;
    fread(., select=c(4,5), col.names=c("REF", "ALT")) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.1.pct
    ## print(after.step.1.pct)

    ## after step 2
    glue("{sample.specific.S15.deeper.path}/alignment.ref.vcf") -> .;
    fread(., header=FALSE, select=c(4,5), col.names=c("REF", "ALT")) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.2.pct
    ## print(after.step.2.pct)

    ## after step 3
    glue("{sample.specific.S15.deeper.path}/alignment.rem.vcf") -> .;
    fread(., header=FALSE, select=c(4,5), col.names=c("REF", "ALT")) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.3.pct
    ## print(after.step.3.pct)

    ## after step 4, Alu subset
    glue("{sample.specific.S15.deeper.path}/alignment.Alu.vcf") -> .;
    fread(., header=FALSE, select=c(4,5), col.names=c("REF", "ALT")) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.4.Alu.pct
    ## print(after.step.4.Alu.pct)

    ## after step 4, other subset
    glue("{sample.specific.S15.deeper.path}/alignment.others.vcf") -> .;
    fread(., header=FALSE, select=c(4,5), col.names=c("REF", "ALT")) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.4.other.pct
    ## print(after.step.4.other.pct)

    ## after step 5
    glue("{sample.specific.S15.deeper.path}/alignment.others.sim.vcf") -> .;
    fread(., header=FALSE, select=c(4,5), col.names=c("REF", "ALT")) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.5.pct
    ## print(after.step.5.pct)

    ## after step 6
    glue("{sample.specific.S15.deeper.path}/alignment.others.rmSJandHomo.txt") -> .;
    fread(., header=FALSE, select=c(4,5), col.names=c("REF", "ALT")) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.6.pct
    ## print(after.step.6.pct)

    ## after step 7
    glue("{sample.specific.S15.deeper.path}/alignment.others.blat.vcf") -> .;
    fread(., header=FALSE, select=c(4,5), col.names=c("REF", "ALT")) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.7.pct
    ## print(after.step.7.pct)

    ## after step 8, non-Alu:
    glue("{sample.specific.S15.deeper.path}/alignment.RepNOTAlu.vcf") -> .;
    fread(., header=FALSE, select=c(4,5), col.names=c("REF", "ALT")) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.8.not.Alu.pct
    ## print(after.step.8.not.Alu.pct)

    ## after step 8, not-repeat:
    glue("{sample.specific.S15.deeper.path}/alignment.nonRep.vcf") -> .;
    fread(., header=FALSE, select=c(4,5), col.names=c("REF", "ALT")) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.8.non.Rep.pct
    ## print(after.step.8.non.Rep.pct)

    ## after step 9
    glue("cat {sample.specific.S15.deeper.path}/alignment.all.real.rich.vcf | grep -v '^##' ") -> .;
    fread(cmd=., header=TRUE) -> .
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(paste(sep="", REF, ">", ALT)) %in% c("A>G", "T>C"))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.9.pct
    ## print(after.step.9.pct)

    ## after step 10
    "result/S51_3__filter_for_variants_with_enough_read_support/210215-sixth-dataset/merged.long.disjoint.with.population.without.potential.polymorphism.dt.txt.gz" -> .;
    glue("zcat {.} | grep ',{sample.name},'") -> .;
    fread(cmd=., header=FALSE, select=1, col.names="ID") -> .;
    .[, type:=sub(pattern=".*_.*_(.*)_(.*)", replacement="\\1>\\2", x=ID)] -> .;
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(type %in% c("A>G", "T>C")))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.10.pct
    ## print(after.step.10.pct)

    ## after step 11
    "result/S51_3__filter_for_variants_with_enough_read_support/210215-sixth-dataset/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt.txt.gz" -> .;
    glue("zcat {.} | grep ',{sample.name},'") -> .;
    fread(cmd=., header=FALSE, select=1, col.names="ID") -> .;
    .[, type:=sub(pattern=".*_.*_(.*)_(.*)", replacement="\\1>\\2", x=ID)] -> .;
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(type %in% c("A>G", "T>C")))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.11.pct
    ## print(after.step.11.pct)

    ## after step 12
    "result/S51_4__filter_for_variants_with_enough_sample_support/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt.txt.gz" -> .;
    glue("zcat {.} | grep '^{sample.name},'") -> .;
    fread(cmd=., header=FALSE, select=2, col.names="ID") -> .;
    .[, type:=sub(pattern=".*_.*_(.*)_(.*)", replacement="\\1>\\2", x=ID)] -> .;
    .[, list(count=.N), list(is.A.to.G.or.T.to.C=(type %in% c("A>G", "T>C")))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or.T.to.C == TRUE, pct] -> after.step.12.pct
    ## print(after.step.12.pct)


    ## after step 13
    "result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz" -> .;
    glue("zcat {.} | grep ',{sample.name},'") -> .;
    fread(cmd=., header=FALSE, select=c(1, 25), col.names=c("ID", "event.summary")) -> .;
    .[, list(count=.N), list(is.A.to.G.or._A.to.G.and.T.to.C_=(event.summary %in% c("A>G", "A>G;T>C")))] -> .;
    .[, pct:=count/sum(count)*100] -> .;
    .[is.A.to.G.or._A.to.G.and.T.to.C_ == TRUE, pct] -> during.step.13.pct
    ## print(during.step.13.pct)
    
}

print(data.table(
    GATK.pct=raw.pct,
    after.step.1.pct=after.step.1.pct,
    after.step.2.pct=after.step.2.pct,
    after.step.3.pct=after.step.3.pct,
    after.step.4.Alu.pct=after.step.4.Alu.pct,
    after.step.4.other.pct=after.step.4.other.pct,
    after.step.5.pct=after.step.5.pct,
    after.step.6.pct=after.step.6.pct,
    after.step.7.pct=after.step.7.pct,
    after.step.8.not.Alu.pct=after.step.8.not.Alu.pct,
    after.step.8.non.Rep.pct=after.step.8.non.Rep.pct,
    after.step.9.pct=after.step.9.pct,
    after.step.10.pct=after.step.10.pct,
    after.step.11.pct=after.step.11.pct,
    after.step.12.pct=after.step.12.pct,
    during.step.13.pct=during.step.13.pct
))


