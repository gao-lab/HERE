library("data.table")
library("stringr")
library("foreach")
library("iterators")

input.gencode.3utr.and.edit.CJ.dt.csv.gz.filename <- snakemake@input[["gencode_3utr_and_edit_CJ_dt_csv_gz_filename"]]
input.gencode.3utr.maf.table.with.species.id.for.a.single.chromosome.filename <- snakemake@input[["gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename"]]
output.edited.gencode.3utr.maf.table.with.species.id.for.a.single.chromosome.filename <- snakemake@output[["edited_gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename"]]

gencode.3utr.and.edit.CJ.dt <- fread(input.gencode.3utr.and.edit.CJ.dt.csv.gz.filename)
raw.3UTR.alignment.dt <- fread(input.gencode.3utr.maf.table.with.species.id.for.a.single.chromosome.filename, header=FALSE, col.names=c("transcript.id", "id", "alignment"))


gencode.3utr.and.edit.edited.block.only.with.human.alignment.dt <- merge(
    x=gencode.3utr.and.edit.CJ.dt[edit.is.in.UTR.block == TRUE],
    y=raw.3UTR.alignment.dt[id==9606],
    by.x=c("transcript_id"),
    by.y=c("transcript.id"),
    all.x=FALSE, all.y=FALSE)

gencode.3utr.and.edit.edited.block.only.with.human.alignment.dt[, edited.alignment:=lapply(alignment, FUN=function(temp.alignment){
        temp.alignment.byte.vector <- charToRaw(temp.alignment)
        temp.to.edit.base <- rawToChar(temp.alignment.byte.vector[edit.rel.POS.wrt.3UTR])
        temp.edited.base.byte <- charToRaw("_")
        if (temp.to.edit.base %in% c("A", "a") == FALSE){x
            warning("Base is not A/a: ", temp.to.edit.base, ", for transcript ", transcript_id, " @ edited at ", edit.POS, " / (relative) ", edit.rel.POS.wrt.3UTR)
        } else {
            temp.edited.base.byte <- c("A"="G", "a"="g")[temp.to.edit.base] %>% charToRaw
        }
        temp.edited.alignment <- temp.alignment.byte.vector %>% {.[edit.rel.POS.wrt.3UTR] <- temp.edited.base.byte; .} %>% rawToChar
        temp.edited.alignment
    }), list(transcript_id, edit.POS)]


## skip edit-free chromosomes (here it is chrY)
if (nrow(gencode.3utr.and.edit.edited.block.only.with.human.alignment.dt) == 0) {
    warning("No 3'-UTR edited in this file: ", input.gencode.3utr.maf.table.with.species.id.for.a.single.chromosome.filename)
    writeLines("NA", output.edited.gencode.3utr.maf.table.with.species.id.for.a.single.chromosome.filename)
    q(save="no")
}



## stack all alignments for each editing
## 211108: decided to use foreach instead of merge to make the syntax easier to understand

gencode.3utr.and.edit.edited.block.only.with.human.alignment.dt[, transcript_and_edit_info:=paste(sep="", transcript_id, "@", edit.POS, "@", edit.rel.POS.wrt.3UTR)]

edited.3UTR.alignment.dt <- foreach(temp.row.dt=iter(gencode.3utr.and.edit.edited.block.only.with.human.alignment.dt, by="row")) %do% {
    temp.raw.3UTR.alignment.subset.dt <- raw.3UTR.alignment.dt[transcript.id == temp.row.dt[1, transcript_id]]
    temp.raw.3UTR.alignment.subset.nonhuman.dt <- temp.raw.3UTR.alignment.subset.dt[id != 9606]
    temp.edited.3UTR.alignment.dt <- data.table(transcript_and_edit_info=temp.row.dt[1, transcript_and_edit_info], rbindlist(list(
        temp.row.dt[, list(id=9606, final.alignment=edited.alignment)],
        temp.raw.3UTR.alignment.subset.nonhuman.dt[, list(id, final.alignment=alignment)]
        )))
    temp.edited.3UTR.alignment.dt
} %>% rbindlist

fwrite(edited.3UTR.alignment.dt, output.edited.gencode.3utr.maf.table.with.species.id.for.a.single.chromosome.filename, sep="\t", col.names=FALSE, quote=FALSE)
