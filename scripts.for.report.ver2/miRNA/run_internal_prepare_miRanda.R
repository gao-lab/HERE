library("data.table")
library("magrittr")
library("foreach")
library("iterators")
library("glue")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/"
dir.create(output.directory, recursive=TRUE)


if (FALSE){
    
    input.original.miRanda.output.gz.filename <- "result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step02__run_miRanda_per_chromosome/32/gencode.3utr.chr22.miRanda.output.test.small.gz"
    input.edited.miRanda.output.gz.filename <- "result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step04__run_miRanda_per_edited_chromosome/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/32/edited.gencode.3utr.chr22.miRanda.output.test.small.gz"
    
}


input.original.miRanda.output.gz.filename <- "result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step05__concatenate_miRanda_results_across_all_chromosomes/32/gencode.3utr.all.chromosomes.concatenated.headless.miRanda.output.gz"
input.edited.miRanda.output.gz.filename <- "result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step06__concatenate_edited_miRanda_results_across_all_chromosomes/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/32/edited.gencode.3utr.all.chromosomes.but.chrY.concatenated.headless.miRanda.output.gz"

input.reference.GTF.filename <- "external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf"
gene.transcript.mapping.dt <- fread(cmd=paste(sep="", "cat ", input.reference.GTF.filename, " | grep -v '^#' | grep -P '\ttranscript\t' | cut -f 9 | sed -E -e 's@.*gene_id \"([^\"]+)\"; transcript_id \"([^\"]+)\";.*@\\1\t\\2@' "), header=FALSE, sep="\t", col.names=c("gene.id", "transcript.id"))


## 1. get original miRanda results
original.miRanda.human.dt <- {
    ##
    input.original.miRanda.output.gz.filename ->.;
    fread(., header=FALSE, sep="\t", select=c(1,2,7), col.names=c("miRNA.family.ID", "transcript.id.3UTR", "UTR.start.and.end")) -> .;
    ## formalize columns
    .[, list(
        transcript.id=sub(pattern="__3UTR$", replacement="", x=transcript.id.3UTR),
        miRNA.family.ID=miRNA.family.ID,
        UTR.start=as.integer(sub(pattern="([0-9]+) ([0-9]+)", replacement="\\1", x=UTR.start.and.end)),
        UTR.end=as.integer(sub(pattern="([0-9]+) ([0-9]+)", replacement="\\2", x=UTR.start.and.end))
    )] -> .;
    setkey(., "transcript.id")
    ##
    .
}
fwrite(original.miRanda.human.dt, glue("{output.directory}/original.miRanda.human.dt.gz"))

## 2. get edited miRanda results
edited.miRanda.human.dt <- {
    ##
    input.edited.miRanda.output.gz.filename -> .;
    fread(., header=FALSE, sep="\t", select=c(1,2,7), col.names=c("miRNA.family.ID", "transcript.and.edit.POS.id.3UTR", "UTR.start.and.end")) -> .;    
    ##= extract columns and formalize column names
    .[, list(
        transcript.id=sub(pattern="^([^@]+)@.*", replacement="\\1", x=transcript.and.edit.POS.id.3UTR),
        edit.POS=sub(pattern="^[^@]+@([^@]+)@.*", replacement="\\1", x=transcript.and.edit.POS.id.3UTR) %>% as.integer,
        edit.rel.POS.wrt.3UTR=sub(pattern="^[^@]+@[^@]+@([^@]+)__3UTR$", replacement="\\1", x=transcript.and.edit.POS.id.3UTR) %>% as.integer,
        miRNA.family.ID=miRNA.family.ID,
        UTR.start=as.integer(sub(pattern="([0-9]+) ([0-9]+)", replacement="\\1", x=UTR.start.and.end)),
        UTR.end=as.integer(sub(pattern="([0-9]+) ([0-9]+)", replacement="\\2", x=UTR.start.and.end))
    )] -> .;
    ##= setkey
    setkey(., "transcript.id", "edit.POS")
    ##
    .
}
fwrite(edited.miRanda.human.dt, glue("{output.directory}/edited.miRanda.human.dt.gz"))



## 3. merge the two to get a list of gain/lost miRNA binding sites per transcript x edit POS
#### 3.1. base table (transcript x edit)
edited.miRanda.human.compared.with.original.dt <- {
    ##
    copy(edited.miRanda.human.dt) -> .;
    ##c get all combinations of transcript X edit POS
    temp.all.pairs.of.transcript.v.s.edit.POS.dt <- unique(.[, list(transcript.id, edit.POS)])
    ##
    print(glue("Total number of records to process: {nrow(temp.all.pairs.of.transcript.v.s.edit.POS.dt)}"))
    foreach(temp.row.dt=iter(temp.all.pairs.of.transcript.v.s.edit.POS.dt[, index:=.I], by="row")) %do% {
        temp.index <- temp.row.dt[1, index]
        if (temp.index %% 1000 == 0){
            cat(date(), " : processing index ", temp.index, "\n")
        }
        temp.combination.dt <- temp.row.dt[, list(transcript.id, edit.POS)]
        temp.original.miRanda.human.subset.dt <- original.miRanda.human.dt[temp.combination.dt[1, transcript.id]]
        temp.edited.miRanda.human.subset.dt <- edited.miRanda.human.dt[temp.combination.dt]
        ## merge all columns but transcript and edit info first
        temp.merged.dt <- merge(
            x=temp.original.miRanda.human.subset.dt[, list(miRNA.family.site.count.from.original=.N), list(miRNA.family.ID)],
            y=temp.edited.miRanda.human.subset.dt[, list(miRNA.family.site.count.from.edited=.N), list(miRNA.family.ID)],
            by=c("miRNA.family.ID"), all.x=TRUE, all.y=TRUE)
        temp.merged.dt[is.na(miRNA.family.site.count.from.original)==TRUE, miRNA.family.site.count.from.original:=0]
        temp.merged.dt[is.na(miRNA.family.site.count.from.edited)==TRUE, miRNA.family.site.count.from.edited:=0]
        ## add transcript and edit info afterwards here so that we know which edit lost which miRNA binding sites
        data.table(temp.edited.miRanda.human.subset.dt[1, list(transcript.id, edit.POS, edit.rel.POS.wrt.3UTR)], temp.merged.dt)
    } -> .; rbindlist(., use.names=TRUE) ->.;
    ##
    .
}

fwrite(edited.miRanda.human.compared.with.original.dt, glue("{output.directory}/edited.miRanda.human.compared.with.original.dt.csv.gz"))

#### 3.2. annotated table (transcript x edit, with annotation)
edited.miRanda.human.compared.with.original.annotated.dt <- {
    ##
    copy(edited.miRanda.human.compared.with.original.dt) -> .;
    ##= add miRNA site affected type    
    .[miRNA.family.site.count.from.original == miRNA.family.site.count.from.edited, miRNA.site.affected.type:="unchanged"]
    .[miRNA.family.site.count.from.original < miRNA.family.site.count.from.edited, miRNA.site.affected.type:="gained"]
    .[miRNA.family.site.count.from.original > miRNA.family.site.count.from.edited, miRNA.site.affected.type:="lost"]    
    ##= add miRNA family ID affected type
    .[miRNA.family.site.count.from.original > 0 & miRNA.family.site.count.from.edited > 0, miRNA.family.ID.affected.type:="unchanged"]
    .[miRNA.family.site.count.from.original == 0 & miRNA.family.site.count.from.edited > 0, miRNA.family.ID.affected.type:="gained"]
    .[miRNA.family.site.count.from.original > 0 & miRNA.family.site.count.from.edited == 0, miRNA.family.ID.affected.type:="lost"]
    ##
    .
}

fwrite(edited.miRanda.human.compared.with.original.annotated.dt, glue("{output.directory}/edited.miRanda.human.compared.with.original.annotated.dt.csv.gz"))

#### 3.3. summarized table (transcript x edit, with net effect of each edit on each transcript)
##c compute difference of count of miRNA sites and families PER transcript X edit 
edited.miRanda.human.compared.with.original.annotated.summary.dt <- edited.miRanda.human.compared.with.original.annotated.dt[
, list(
      difference.of.count.of.miRNA.sites=sum(miRNA.family.site.count.from.edited - miRNA.family.site.count.from.original),
      difference.of.count.of.miRNA.families=sum(miRNA.family.ID.affected.type=='gained') - sum(miRNA.family.ID.affected.type=='lost')
  ), list(transcript.id, edit.POS)]

fwrite(edited.miRanda.human.compared.with.original.annotated.summary.dt, glue("{output.directory}/edited.miRanda.human.compared.with.original.annotated.summary.dt.csv.gz"))

#### 3.4. summarized table at gene level (gene x edit, with net effect of each edit on each gene)
edited.miRanda.human.compared.with.original.annotated.summary.gene.and.edit.level.dt <- {
    copy(edited.miRanda.human.compared.with.original.annotated.summary.dt) ->.;
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

fwrite(edited.miRanda.human.compared.with.original.annotated.summary.gene.and.edit.level.dt, glue("{output.directory}/edited.miRanda.human.compared.with.original.annotated.summary.gene.and.edit.level.dt.csv.gz"))


