library("glue")
library("data.table")
library("foreach")
library("iterators")
library("doMC")
library("magrittr")
library("readxl")
source("./scripts/common/ggpubr.A4.R")



output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/lost.edits.coverage.check/"
dir.create(output.directory, recursive=TRUE)

## TTF1 and ELP3 part
{

    ## 1.1. select all normal zygotes and all GSE133854 zygotes samples
    {
        fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt") -> .;
        ## keep normal early stage samples
        .[stage == 'zygote'][(is.normal == TRUE) | (gse=="GSE133854")] -> .;
        .
    } -> ..normal.or.GSE133854.zygote.samples.info.dt

    ## 1.2. link bams
    {
        copy(..normal.or.GSE133854.zygote.samples.info.dt) -> .;
        .[gse=="GSE71318", original.bam.path:=paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-125-125/200919-GSE71318-full48-125-125/", gsm, "/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/120/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam")] -> .;
        .[gse=="GSE133854" & srr.mean.avgspotlen==180, original.bam.path:=paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-90-90/200924-GSE133854-all296-90-90/", gsm, "/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/85/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam")] -> .;
        .[gse=="GSE133854" & srr.mean.avgspotlen==300, original.bam.path:=paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-150-150/200924-GSE133854-all296-150-150/", gsm, "/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/145/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam")] -> .;
        .[gse=="GSE44183", original.bam.path:=paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-90-90/201217-GSE44183-earlyhumanlong21-90-90/", gsm, "/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/85/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam")] -> .;
        .[gse=="GSE36552", original.bam.path:=paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/single-100/201104-GSE36552-full124-100/", gsm, "/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/95/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam")] -> .;
        ##.[, original.bam.path:=glue("../Project.collaborated/Project.embryo.I.201220/{original.bam.path}")] -> .;
        .[, original.bai.path:=sub(pattern="recal\\.bam$", replacement="recal\\.bai", x=original.bam.path)] -> .;
        
        ## build link
        foreach(temp.row.dt=iter(., by="row")) %do% {
            ..temp.gsm.directory <- glue("{output.directory}/{temp.row.dt[1, gsm]}")
            dir.create(..temp.gsm.directory, recursive=TRUE)
            system(glue("rm {..temp.gsm.directory}/temp.recal.bam"))
            system(glue("ln -s -r {temp.row.dt[1, original.bam.path]} {..temp.gsm.directory}/temp.recal.bam"))
            system(glue("rm {..temp.gsm.directory}/temp.recal.bai"))
            system(glue("ln -s -r {temp.row.dt[1, original.bai.path]} {..temp.gsm.directory}/temp.recal.bai"))            
        }
        
    } 

    ## 1.3. subset linked bams
    ## chr9, 132375956
    {

        ## prepare for the edit site
        ..TTF1.edit.bed.filename <- glue("{output.directory}/TTF1.edit.bed")
        fwrite(data.table(chrom="chr9", chromStart=132375956 - 1, chromEnd=132375956 - 1 + 1), ..TTF1.edit.bed.filename, row.names=FALSE, col.names=FALSE, sep="\t")

        ## subset linked bams using the edit site
        foreach(temp.row.dt=iter(..normal.or.GSE133854.zygote.samples.info.dt, by="row")) %do% {
            cat(glue("{date()} : processing sample {temp.row.dt[1, gsm]} "), "\n")
            ..temp.gsm.directory <- glue("{output.directory}/{temp.row.dt[1, gsm]}")
            system(glue("samtools view -@ 8 -L {..TTF1.edit.bed.filename} -o {..temp.gsm.directory}/temp.recal.chr9.132375956.subset.bam {..temp.gsm.directory}/temp.recal.bam "))
            system(glue("samtools index {..temp.gsm.directory}/temp.recal.chr9.132375956.subset.bam"))
        }

    }


    ## chr8, 28190741
    {

        ## prepare for the edit site
        ..ELP3.edit.bed.filename <- glue("{output.directory}/ELP3.edit.bed")
        fwrite(data.table(chrom="chr8", chromStart=28190741 - 1, chromEnd=28190741 - 1 + 1), ..ELP3.edit.bed.filename, row.names=FALSE, col.names=FALSE, sep="\t")

        ## subset linked bams using the edit site
        foreach(temp.row.dt=iter(..normal.or.GSE133854.zygote.samples.info.dt, by="row")) %do% {
            cat(glue("{date()} : processing sample {temp.row.dt[1, gsm]} "), "\n")
            ..temp.gsm.directory <- glue("{output.directory}/{temp.row.dt[1, gsm]}")
            system(glue("samtools view -@ 8 -L {..ELP3.edit.bed.filename} -o {..temp.gsm.directory}/temp.recal.chr8.28190741.subset.bam {..temp.gsm.directory}/temp.recal.bam "))
            system(glue("samtools index {..temp.gsm.directory}/temp.recal.chr8.28190741.subset.bam"))
        }

    }

    ## 1.4. rename result bams (such that they need no renaming in IGV)
    ## 1.1. get editing level for each edited events on valid genes in normal early stages
    {
        
        ## read edits in all samples
        fread(
            "result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz"
        ) -> ..total.dt;

        ## TTF1 part
        ## pick normal or GSE133854 zygotes
        ..total.dt[stage == 'zygote'][(is.normal == TRUE) | (gse == "GSE133854")] -> .;
        ## pick the edit chr9:132375956
        .[CHROM=='chr9' & POS==132375956] -> .;
        ## merge the edit info with ..normal.or.GSE133854.zygote.samples.info.dt
        merge(x=..normal.or.GSE133854.zygote.samples.info.dt, y=.[, list(SAMPLE, AC, AN, AF)],
              by.x="gsm", by.y="SAMPLE",
              all.x=TRUE, all.y=FALSE) -> .;
        .[, detection.type:="detected"][is.na(AC)==TRUE, detection.type:="undetected"] -> .;
        .[, status:="normal-other"][(gse=="GSE133854") & ((disease %in% c("androgenetic", "pathenogenetic")) == FALSE), status:="BI"][disease=="androgenetic", status:="AG"][disease=="parthenogenetic", status:="PG"] -> .;
        .[, new.bam.basename:=paste(sep="", detection.type, "_", gse, "_", gsm, "_", status, "_", AF, "_", AC, "_", AN, ".bam")] -> .;
        .[, new.bai.basename:=sub(pattern="\\.bam$", replacement=".bam.bai", x=new.bam.basename)] -> .;
        ## cp the subset bams to new names
        foreach(temp.row.dt=iter(., by="row")) %do% {
            cat(glue("{date()} : processing sample {temp.row.dt[1, gsm]} "), "\n")
            ..temp.gsm.directory <- glue("{output.directory}/{temp.row.dt[1, gsm]}")
            system(glue("cp {..temp.gsm.directory}/temp.recal.chr9.132375956.subset.bam {..temp.gsm.directory}/{temp.row.dt[1, new.bam.basename]}"))
            system(glue("cp {..temp.gsm.directory}/temp.recal.chr9.132375956.subset.bam.bai {..temp.gsm.directory}/{temp.row.dt[1, new.bai.basename]}"))            
        }


        ## ELP3 part
        ## pick normal or GSE133854 zygotes
        ..total.dt[stage == 'zygote'][(is.normal == TRUE) | (gse == "GSE133854")] -> .;
        ## pick the edit chr8:28190741
        .[CHROM=='chr8' & POS==28190741] -> .;
        ## merge the edit info with ..normal.or.GSE133854.zygote.samples.info.dt
        merge(x=..normal.or.GSE133854.zygote.samples.info.dt, y=.[, list(SAMPLE, AC, AN, AF)],
              by.x="gsm", by.y="SAMPLE",
              all.x=TRUE, all.y=FALSE) -> .;
        .[, detection.type:="detected"][is.na(AC)==TRUE, detection.type:="undetected"] -> .;
        .[, status:="normal-other"][(gse=="GSE133854") & ((disease %in% c("androgenetic", "pathenogenetic")) == FALSE), status:="BI"][disease=="androgenetic", status:="AG"][disease=="parthenogenetic", status:="PG"] -> .;
        .[, new.bam.basename:=paste(sep="", "ELP3_chr8_28190741__", detection.type, "_", gse, "_", gsm, "_", status, "_", AF, "_", AC, "_", AN, ".bam")] -> .;
        .[, new.bai.basename:=sub(pattern="\\.bam$", replacement=".bam.bai", x=new.bam.basename)] -> .;
        ## cp the subset bams to new names
        foreach(temp.row.dt=iter(., by="row")) %do% {
            cat(glue("{date()} : processing sample {temp.row.dt[1, gsm]} "), "\n")
            ..temp.gsm.directory <- glue("{output.directory}/{temp.row.dt[1, gsm]}")
            system(glue("cp {..temp.gsm.directory}/temp.recal.chr8.28190741.subset.bam {..temp.gsm.directory}/{temp.row.dt[1, new.bam.basename]}"))
            system(glue("cp {..temp.gsm.directory}/temp.recal.chr8.28190741.subset.bam.bai {..temp.gsm.directory}/{temp.row.dt[1, new.bai.basename]}"))            
        }

    }
    
}



## 107 edits part
{

    ## 2.1. select all GSE133854 abnormal and GSE95477 elder mother samples
    {
        fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt") -> .;
        ## keep GSE133854 abnormal and GSE95477 elder mother samples
        .[gsm %in% union(.[gse=="GSE95477"][as.integer(maternal.age) > 35][, gsm], .[gse=="GSE133854"][is.normal==FALSE][, gsm])] -> .;
        .
    } -> ..GSE95477.elder.mother.or.GSE133854.abnormal.samples.info.dt

    ## 2.2. link bams
    {
        
        copy(..GSE95477.elder.mother.or.GSE133854.abnormal.samples.info.dt) -> .;
        .[gse=="GSE95477", original.bam.path:=paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-100-100/201101-GSE95477-full20-100-100/", gsm, "/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/95/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam")] -> .;
        .[gse=="GSE133854" & srr.mean.avgspotlen==180, original.bam.path:=paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-90-90/200924-GSE133854-all296-90-90/", gsm, "/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/85/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam")] -> .;
        .[gse=="GSE133854" & srr.mean.avgspotlen==300, original.bam.path:=paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-150-150/200924-GSE133854-all296-150-150/", gsm, "/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/145/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam")] -> .;
        ##.[, original.bam.path:=glue("../Project.collaborated/Project.embryo.I.201220/{original.bam.path}")] -> .;
        .[, original.bai.path:=sub(pattern="recal\\.bam$", replacement="recal\\.bai", x=original.bam.path)] -> .;
        
        ## build link
        foreach(temp.row.dt=iter(., by="row")) %do% {
            ..temp.gsm.directory <- glue("{output.directory}/107edits/{temp.row.dt[1, gsm]}")
            dir.create(..temp.gsm.directory, recursive=TRUE)
            system(glue("rm {..temp.gsm.directory}/temp.recal.bam"))
            system(glue("ln -s -r {temp.row.dt[1, original.bam.path]} {..temp.gsm.directory}/temp.recal.bam"))
            system(glue("rm {..temp.gsm.directory}/temp.recal.bai"))
            system(glue("ln -s -r {temp.row.dt[1, original.bai.path]} {..temp.gsm.directory}/temp.recal.bai"))            
        }
        
    } 

    ## 2.3. generate bed file to compute coverage
    {
        "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/disease.or.old.mother.embryo.lost.RE.dt.csv.gz" -> .;
        fread(.) -> .;
        ## generate bed file
        fwrite(unique(.[, list(chrom=CHROM, chromStart=POS-1, chromEnd=POS-1+1)]), glue("{output.directory}/107.edits.bed"), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
    }

    
    ## 2.4. compute the sequencing coverage using the bed
    {

        registerDoMC(cores=10)
        foreach(temp.row.dt=iter(..GSE95477.elder.mother.or.GSE133854.abnormal.samples.info.dt, by="row")) %dopar% {
            cat(glue("{date()} : processing sample {temp.row.dt[1, gsm]} "), "\n")
            ..temp.gsm.directory <- glue("{output.directory}/107edits/{temp.row.dt[1, gsm]}")
            system(glue("samtools mpileup --ignore-overlaps --min-BQ=0 -l {output.directory}/107.edits.bed -o {..temp.gsm.directory}/temp.recal.107editsbed.pileup {..temp.gsm.directory}/temp.recal.bam "))
        }

        rbindlist(foreach(temp.row.dt=iter(..GSE95477.elder.mother.or.GSE133854.abnormal.samples.info.dt, by="row")) %do% {
            ..temp.gsm.directory <- glue("{output.directory}/107edits/{temp.row.dt[1, gsm]}")
            fread(glue("{..temp.gsm.directory}/temp.recal.107editsbed.pileup"), header=FALSE, col.names=c("CHROM", "POS", "ref", "count.of.bases", "read.bases", "base.qualities"))[, gsm:=temp.row.dt[1, gsm]]
        }, use.names=TRUE) -> .;
        .
    } -> ..all.107edits.pileup.dt
    fwrite(..all.107edits.pileup.dt, glue("{output.directory}/107edits/all.107edits.pileup.dt.gz"))

    ## 2.5. generate the per-edit, per-lost-type coverage table
    {
        
        "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/disease.or.old.mother.embryo.lost.RE.dt.csv.gz" -> .;
        fread(.) -> .;
        unique(.[, list(CHROM, POS, stage, completely.lost.in.GSE95477.old.mother.embryos, completely.lost.in.GSE133854.AG.embryos, completely.lost.in.GSE133854.PG.embryos)]) -> .;
        
        rbindlist(list(
            merge(
                x=.[completely.lost.in.GSE95477.old.mother.embryos == TRUE],
                y=..GSE95477.elder.mother.or.GSE133854.abnormal.samples.info.dt[gse=="GSE95477"],
                by="stage",
                all.x=TRUE, all.y=FALSE, allow.cartesian=TRUE
            ),
            merge(
                x=.[completely.lost.in.GSE133854.AG.embryos == TRUE],
                y=..GSE95477.elder.mother.or.GSE133854.abnormal.samples.info.dt[gse=="GSE133854" & disease == "androgenetic"],
                by="stage",
                all.x=TRUE, all.y=FALSE, allow.cartesian=TRUE
            ),
            merge(
                x=.[completely.lost.in.GSE133854.PG.embryos == TRUE],
                y=..GSE95477.elder.mother.or.GSE133854.abnormal.samples.info.dt[gse=="GSE133854" & disease == "parthenogenetic"],
                by="stage",
                all.x=TRUE, all.y=FALSE, allow.cartesian=TRUE
            )
        ), use.names=TRUE) -> ..all.107edits.and.lost.samples.dt

        merge(
            x=..all.107edits.and.lost.samples.dt[, list(CHROM, POS, gse, gsm, stage, disease)],
            y=..all.107edits.pileup.dt[, list(CHROM, POS, gsm, count.of.bases)],
            by=c("CHROM", "POS", "gsm"),
            all.x=TRUE, all.y=FALSE
        )[is.na(count.of.bases) == TRUE, count.of.bases:=0] -> ..all.107edits.and.lost.sample.count.of.bases.dt

        fwrite(..all.107edits.and.lost.sample.count.of.bases.dt, glue("{output.directory}/all.107edits.and.lost.sample.count.of.bases.dt.gz"))
        ..all.107edits.and.lost.sample.count.of.bases.dt
        
    } -> all.107edits.and.lost.sample.count.of.bases.dt

    ## 2.6. plot the distribution
    {
        
        all.107edits.and.lost.sample.count.of.bases.dt -> .;
        ## prettify stages
        temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table
        merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.reverse.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
        ## prettify label
        .[gse=="GSE133854", label:=c("androgenetic"="AG", "parthenogenetic"="PG")[disease]] -> .;
        .[gse=="GSE95477", label:="elder mother"] -> .;
        ##
        . -> ..to.plot.dt
        
        ggplot(..to.plot.dt[gse=="GSE95477"], aes(x=paste(sep="", CHROM, "_", POS), y=gsm, fill=cut(count.of.bases, breaks=c(-Inf, 0, 9, 99, Inf), labels=c("0", "[1, 9]", "[10, 99]", ">=100")))) ->.;
        . + geom_tile() ->.;
        . + coord_flip() ->.;
        ## facet
        . + facet_wrap(~stage.description.reverse.ordered + label, scales="free", nrow=1) ->.;
        ## add theme
        . + theme_pubr() ->.;
        . + theme(legend.position="top", axis.text.x=element_text(angle=45, hjust=1)) ->.;
        . + labs(x="REE", y="Sample ID", fill="Coverage") ->.;
        . + scale_fill_manual(values=c("#61E2E8", "#77CBFF", "#756BFF")) -> .;
        ## . + scale_y_continuous(breaks=c(0, 15, 30, 45)) -> .;
        ## save image
        ggsave.A4(
            filename=glue("{output.directory}/REE.107edits.sequencing.coverage.in.lost.samples.GSE95477.png"),
            plot=.,
            width.r=0.9, height.r=0.5)


        ggplot(..to.plot.dt[gse=="GSE133854" & stage == 'zygote'], aes(x=paste(sep="", CHROM, "_", POS), y=gsm, fill=cut(count.of.bases, breaks=c(-Inf, 0, 9, 99, Inf), labels=c("0", "[1, 9]", "[10, 99]", ">=100")))) ->.;
        . + geom_tile() ->.;
        . + coord_flip() ->.;
        ## facet
        . + facet_wrap(~stage.description.reverse.ordered + label, scales="free", nrow=1) ->.;
        ## add theme
        . + theme_pubr() ->.;
        . + theme(legend.position="top", axis.text.x=element_text(angle=45, hjust=1)) ->.;
        . + labs(x="REE", y="Sample ID", fill="Coverage") ->.;
        . + scale_fill_manual(values=c("gray50", "#61E2E8", "#77CBFF", "#756BFF")) -> .;
        ## . + scale_y_continuous(breaks=c(0, 15, 30, 45)) -> .;
        ## save image
        ggsave.A4(
            filename=glue("{output.directory}/REE.107edits.sequencing.coverage.in.lost.samples.GSE133854.zygotes.png"),
            plot=.,
            width.r=0.9, height.r=1)

        ggplot(..to.plot.dt[gse=="GSE133854" & stage == '2-cell'], aes(x=paste(sep="", CHROM, "_", POS), y=gsm, fill=cut(count.of.bases, breaks=c(-Inf, 0, 9, 99, Inf), labels=c("0", "[1, 9]", "[10, 99]", ">=100")))) ->.;
        . + geom_tile() ->.;
        . + coord_flip() ->.;
        ## facet
        . + facet_wrap(~stage.description.reverse.ordered + label, scales="free", nrow=1) ->.;
        ## add theme
        . + theme_pubr() ->.;
        . + theme(legend.position="top", axis.text.x=element_text(angle=45, hjust=1)) ->.;
        . + labs(x="REE", y="Sample ID", fill="Coverage") ->.;
        . + scale_fill_manual(values=c("gray50", "#61E2E8", "#77CBFF", "#756BFF")) -> .;
        ## . + scale_y_continuous(breaks=c(0, 15, 30, 45)) -> .;
        ## save image
        ggsave.A4(
            filename=glue("{output.directory}/REE.107edits.sequencing.coverage.in.lost.samples.GSE133854.2cells.png"),
            plot=.,
            width.r=0.9, height.r=0.25)
        
    }
    
    
}
