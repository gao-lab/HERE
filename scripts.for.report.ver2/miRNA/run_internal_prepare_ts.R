library("data.table")
library("magrittr")
library("foreach")
library("iterators")
source("./scripts/common/ggpubr.A4.R")

dir.create("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/", recursive=TRUE)

input.original.targetscan.output.gz.filename <- "result/A02_8__get_editing_effect_on_miRNA_binding_sites/step09__concatenate_TargetScan_results_across_all_chromosomes/32/gencode.3utr.all.chromosomes.concatenated.headless.TargetScan.output.gz"
input.edited.targetscan.output.gz.filename <- "result/A02_8__get_editing_effect_on_miRNA_binding_sites/step10__concatenate_edited_TargetScan_results_across_all_chromosomes/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/32/edited.gencode.3utr.all.chromosomes.but.chrY.concatenated.headless.TargetScan.output.gz"

input.reference.GTF.filename <- "external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf"
gene.transcript.mapping.dt <- fread(cmd=paste(sep="", "cat ", input.reference.GTF.filename, " | grep -v '^#' | grep -P '\ttranscript\t' | cut -f 9 | sed -E -e 's@.*gene_id \"([^\"]+)\"; transcript_id \"([^\"]+)\";.*@\\1\t\\2@' "), header=FALSE, sep="\t", col.names=c("gene.id", "transcript.id"))


## 1. get original TargetScan results
original.ts.human.dt <- {
    input.original.targetscan.output.gz.filename ->.;
    fread(., header=FALSE, col.names=c("a_Gene_ID", "miRNA_family_ID", "species_ID", "MSA_start", "MSA_end", "UTR_start", "UTR_end", "Group_num", "Site_type", "miRNA in this species", "Group_type", "Species_in_this_group", "Species_in_this_group_with_this_site_type", "ORF_overlap")) -> .;
    ##s pick human miRNAs only
    .[species_ID==9606 & `miRNA in this species` == 'x'] ->.;
    ##= formalize column names
    .[, list(
        transcript.id=a_Gene_ID,
        miRNA.family.ID=miRNA_family_ID,
        MSA.start=MSA_start,
        MSA.end=MSA_end,
        UTR.start=UTR_start,
        UTR.end=UTR_end)] -> .;
    ##= setkey
    setkey(., "transcript.id")
}

## 2. get edited TargetScan results
edited.ts.human.dt <- {
    input.edited.targetscan.output.gz.filename -> .;
    fread(., header=FALSE, col.names=c("a_Gene_ID", "miRNA_family_ID", "species_ID", "MSA_start", "MSA_end", "UTR_start", "UTR_end", "Group_num", "Site_type", "miRNA in this species", "Group_type", "Species_in_this_group", "Species_in_this_group_with_this_site_type", "ORF_overlap")) ->.;
    ##s pick human miRNAs only
    .[species_ID==9606 & `miRNA in this species` == 'x'] -> .;
    ##= extract columns and formalize column names
    .[, list(
        transcript.id=sub(pattern="^([^@]+)@.*", replacement="\\1", x=a_Gene_ID),
        edit.POS=sub(pattern="^[^@]+@([^@]+)@.*", replacement="\\1", x=a_Gene_ID) %>% as.integer,
        edit.rel.POS.wrt.3UTR=sub(pattern="^[^@]+@[^@]+@([^@]+)$", replacement="\\1", x=a_Gene_ID) %>% as.integer,
        miRNA.family.ID=miRNA_family_ID,
        MSA.start=MSA_start,
        MSA.end=MSA_end,
        UTR.start=UTR_start,
        UTR.end=UTR_end)] ->.;
    ##= setkey
    setkey(., "transcript.id", "edit.POS")
}


## 3. merge the two to get a list of gain/lost miRNA binding sites per transcript x edit POS
#### 3.1. base table (transcript x edit)
edited.ts.human.compared.with.original.dt <- {
    copy(edited.ts.human.dt) -> .;
    ##c get all combinations of transcript X edit POS
    temp.all.pairs.of.transcript.v.s.edit.POS.dt <- unique(.[, list(transcript.id, edit.POS)])
    foreach(temp.row.dt=iter(temp.all.pairs.of.transcript.v.s.edit.POS.dt[, index:=.I], by="row")) %do% {
        temp.index <- temp.row.dt[1, index]
        if (temp.index %% 1000 == 0){
            cat(date(), " : processing index ", temp.index, "\n")
        }
        temp.combination.dt <- temp.row.dt[, list(transcript.id, edit.POS)]
        temp.original.ts.human.subset.dt <- original.ts.human.dt[temp.combination.dt[1, transcript.id]]
        temp.edited.ts.human.subset.dt <- edited.ts.human.dt[temp.combination.dt]
        ## merge all columns but transcript and edit info first
        temp.merged.dt <- merge(
            x=temp.original.ts.human.subset.dt[, list(miRNA.family.site.count.from.original=.N), list(miRNA.family.ID)],
            y=temp.edited.ts.human.subset.dt[, list(miRNA.family.site.count.from.edited=.N), list(miRNA.family.ID)],
            by=c("miRNA.family.ID"), all.x=TRUE, all.y=TRUE)
        temp.merged.dt[is.na(miRNA.family.site.count.from.original)==TRUE, miRNA.family.site.count.from.original:=0]
        temp.merged.dt[is.na(miRNA.family.site.count.from.edited)==TRUE, miRNA.family.site.count.from.edited:=0]
        ## add transcript and edit info afterwards here so that we know which edit lost which miRNA binding sites
        data.table(temp.edited.ts.human.subset.dt[1, list(transcript.id, edit.POS, edit.rel.POS.wrt.3UTR)], temp.merged.dt)
    } -> .; rbindlist(., use.names=TRUE) ->.;
}

fwrite(edited.ts.human.compared.with.original.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/edited.ts.human.compared.with.original.dt.csv.gz")

#### 3.2. annotated table (transcript x edit, with annotation)
edited.ts.human.compared.with.original.annotated.dt <- {
    copy(edited.ts.human.compared.with.original.dt) -> .;
    ##= add miRNA site affected type    
    .[miRNA.family.site.count.from.original == miRNA.family.site.count.from.edited, miRNA.site.affected.type:="unchanged"]
    .[miRNA.family.site.count.from.original < miRNA.family.site.count.from.edited, miRNA.site.affected.type:="gained"]
    .[miRNA.family.site.count.from.original > miRNA.family.site.count.from.edited, miRNA.site.affected.type:="lost"]    
    ##= add miRNA family ID affected type
    .[miRNA.family.site.count.from.original > 0 & miRNA.family.site.count.from.edited > 0, miRNA.family.ID.affected.type:="unchanged"]
    .[miRNA.family.site.count.from.original == 0 & miRNA.family.site.count.from.edited > 0, miRNA.family.ID.affected.type:="gained"]
    .[miRNA.family.site.count.from.original > 0 & miRNA.family.site.count.from.edited == 0, miRNA.family.ID.affected.type:="lost"]
}

fwrite(edited.ts.human.compared.with.original.annotated.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/edited.ts.human.compared.with.original.annotated.dt.csv.gz")

#### 3.3. summarized table (transcript x edit, with net effect of each edit on each transcript)
##c compute difference of count of miRNA sites and families PER transcript X edit 
edited.ts.human.compared.with.original.annotated.summary.dt <- edited.ts.human.compared.with.original.annotated.dt[
, list(
      difference.of.count.of.miRNA.sites=sum(miRNA.family.site.count.from.edited - miRNA.family.site.count.from.original),
      difference.of.count.of.miRNA.families=sum(miRNA.family.ID.affected.type=='gained') - sum(miRNA.family.ID.affected.type=='lost')
  ), list(transcript.id, edit.POS)]

fwrite(edited.ts.human.compared.with.original.annotated.summary.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/edited.ts.human.compared.with.original.annotated.summary.dt.csv.gz")

#### 3.4. summarized table at gene level (gene x edit, with net effect of each edit on each gene)
edited.ts.human.compared.with.original.annotated.summary.gene.and.edit.level.dt <- {
    copy(edited.ts.human.compared.with.original.annotated.summary.dt) ->.;
    ##m merge with gene info
    merge(x=., y=gene.transcript.mapping.dt,
          by.x="transcript.id", by.y="transcript.id",
          all.x=TRUE, all.y=FALSE) -> .;
    ##c collapse at the level of gene x edit
    .[, list(
        gains.miRNA.sites=all(difference.of.count.of.miRNA.sites >= 0) & any(difference.of.count.of.miRNA.sites > 0),
        gains.miRNA.families=all(difference.of.count.of.miRNA.families >= 0) & any(difference.of.count.of.miRNA.families > 0),
        loses.miRNA.sites=all(difference.of.count.of.miRNA.sites <= 0) & any(difference.of.count.of.miRNA.sites < 0),
        loses.miRNA.families=all(difference.of.count.of.miRNA.families <= 0) & any(difference.of.count.of.miRNA.families < 0)
    ), list(gene.id, edit.POS)] -> .;
}

fwrite(edited.ts.human.compared.with.original.annotated.summary.gene.and.edit.level.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/edited.ts.human.compared.with.original.annotated.summary.gene.and.edit.level.dt.csv.gz")


