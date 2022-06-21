library("data.table")
library("Biostrings")
library("glue")
library("foreach")
library("iterators")
library("stringi")
library("ggpubr")
library("readxl")
library("ggtext")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.motif/"
dir.create(output.directory, recursive=TRUE)

glue("chr{c(1:22, 'X')}") -> chromosomes.to.check.vector
motif.dt <- rbindlist(list(
    data.table(motif.name="Geissler2016", motif.regex="TAA[CG]{1}TTAT"),
    data.table(motif.name="PRE", motif.regex="TGTA[ACGT]{1}AT[AT]{1}"),
    data.table(motif.name="ARE.9mer", motif.regex="TTATTTATT")
))

{

    ## read unedited fasta into data.table
    {
        glue("./result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step01__get_3UTR_sequences_from_maf_blocks/32/gencode.3utr.{chromosomes.to.check.vector}.fasta") -> .
        do.call(c, lapply(., FUN=readDNAStringSet)) -> .
        data.table(raw.name=names(.), sequence=as.character(.), edit.status="original") -> .;
        .
    } -> ..original.sequence.dt;

    ## read edited fasta into data.table
    {
        glue("./result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step03__get_edited_3UTR_sequences_from_edited_maf_blocks/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/32/edited.gencode.3utr.{chromosomes.to.check.vector}.fasta") -> .
        do.call(c, lapply(., FUN=readDNAStringSet)) -> .
        data.table(raw.name=names(.), sequence=as.character(.), edit.status="edited") -> .;
        .
    } -> ..edited.sequence.dt

    ## combine both
    rbindlist(list(..original.sequence.dt, ..edited.sequence.dt), use.names=TRUE) -> ..combined.sequence.dt;
    
    ## get motif position
    ## ~1min per motif
    {
        rbindlist(foreach(TEMP.ROW.DT=iter(motif.dt, by="row")) %do% {
            print(glue("{date()} : processing motif {TEMP.ROW.DT[1, motif.name]}"))
            stri_locate_all_regex(str=..combined.sequence.dt[, sequence], pattern=TEMP.ROW.DT[1, motif.regex]) -> ..temp.stri.hits.list
            rbindlist(mapply(
                ..combined.sequence.dt[, raw.name], ..combined.sequence.dt[, edit.status], ..temp.stri.hits.list,
                FUN=function(..temp.raw.name, ..temp.edit.status, ..temp.positions.matrix){
                    data.table(raw.name=..temp.raw.name, edit.status=..temp.edit.status, ..temp.positions.matrix)
                },
                SIMPLIFY=FALSE
            ), use.names=TRUE) -> TEMP.RESULT.DT;
            TEMP.RESULT.DT[, motif.name:=TEMP.ROW.DT[1, motif.name]] -> TEMP.RESULT.DT;
            TEMP.RESULT.DT
        }, use.names=TRUE) -> .;
        ##
        ## add annotations
        .[
            edit.status=="original",
            transcript.id:=sub(pattern="__3UTR$", replacement="", x=raw.name)
        ] -> .;
        ##
        .[
            edit.status=="edited",
            `:=`(
                transcript.id=sub(pattern="^([^@]+)@.*", replacement="\\1", x=raw.name),
                edit.POS=as.integer(sub(pattern="^[^@]+@([^@]+)@.*", replacement="\\1", x=raw.name)),
                edit.rel.POS.wrt.3UTR=as.integer(sub(pattern="^[^@]+@[^@]+@([^@]+)__3UTR$", replacement="\\1", x=raw.name))
            )
        ] -> .;
        ##
        setnames(., c("start", "end"), c("motif.start", "motif.end")) -> .;
        ##
        .[
        (edit.status=='edited') & (is.na(motif.start) == FALSE),
        the.motif.site.in.the.current.edit.status.has.ge.1.site.overlapping.the.edit := (
            (edit.rel.POS.wrt.3UTR >= motif.start) &
            (edit.rel.POS.wrt.3UTR <= motif.end)
        )] -> .;
        ##
        .
    } -> ..motif.annotation.dt
    
    fwrite(..motif.annotation.dt, glue("{output.directory}/motif.annotation.dt.gz"))
    ..motif.annotation.dt
} -> motif.annotation.dt

## get the list of gain/lost motif per transcript x edit POS
{

    ## [1] get all (transcript.id, edit.POS) pairs
    {
        ##
        "./result/A02_8__get_editing_effect_on_miRNA_binding_sites/step06__compute_edit_relative_position_on_3UTR/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/32/gencode.3utr.and.edit.CJ.dt.csv.gz" -> .;
        fread(.) -> .;
        ## pick edited only
        .[is.na(edit.POS)==FALSE] -> .;
        ## pick those whose edit is in UTR block only
        .[edit.is.in.UTR.block == TRUE] -> .;
        .[, list(transcript.id=transcript_id, edit.POS, edit.rel.POS.wrt.3UTR)] -> .;
        ## add index
        .[, index:=.I] -> .;
        ##
        .
    } -> ..all.transcript.and.edit.pairs.to.check.dt


    ## [2] base table (transcript x edit)    
    {
        ##
        foreach(TEMP.ROW.DT=iter(..all.transcript.and.edit.pairs.to.check.dt, by="row")) %do% {
            TEMP.INDEX <- TEMP.ROW.DT[1, index]
            if (TEMP.INDEX %% 1000 == 0){
                print(glue("{date()} : processing index {TEMP.INDEX}"))
            }
            ##
            ## get original hit counts
            motif.annotation.dt[
            (edit.status=='original') & (transcript.id == TEMP.ROW.DT[1, transcript.id]) & (is.na(motif.start) == FALSE)
            ][,
              the.motif.site.in.the.original.has.ge.1.site.overlapping.the.edit:=(
                  (TEMP.ROW.DT[1, edit.rel.POS.wrt.3UTR] >= motif.start) &
                  (TEMP.ROW.DT[1, edit.rel.POS.wrt.3UTR] <= motif.end)
              )] -> ..temp.original.hit.dt
            ##
            ..temp.original.hit.dt[
              , list(
                    motif.site.count.from.original=.N,
                    from.original=TRUE,
                    the.motif.in.the.original.has.ge.1.site.overlapping.the.edit=any(the.motif.site.in.the.original.has.ge.1.site.overlapping.the.edit)),
                list(motif.name)
            ] -> ..temp.original.hit.count.dt
            ##
            ## get edited hit counts
            motif.annotation.dt[
            (edit.status=='edited') & (transcript.id == TEMP.ROW.DT[1, transcript.id]) & (edit.POS == TEMP.ROW.DT[1, edit.POS]) & (is.na(motif.start)==FALSE)
            ][,
                list(
                    motif.site.count.from.edited=.N,
                    from.edited=TRUE,
                    the.motif.in.the.edited.has.ge.1.site.overlapping.the.edit=any(the.motif.site.in.the.current.edit.status.has.ge.1.site.overlapping.the.edit)),
                list(motif.name)
            ] -> ..temp.edited.hit.count.dt
            ##
            ## merge all columns but transcript and edit info first
            merge(
                x=..temp.original.hit.count.dt, y=..temp.edited.hit.count.dt,
                by=c("motif.name"),
                all.x=TRUE, all.y=TRUE
            ) -> ..temp.merged.dt
            ##
            ## fill the zeroes and FALSE's
            ..temp.merged.dt[is.na(from.original)==TRUE, `:=`(motif.site.count.from.original=0, from.original=FALSE)] -> ..temp.merged.dt;
            ..temp.merged.dt[is.na(from.edited)==TRUE, `:=`(motif.site.count.from.edited=0, from.edited=FALSE)] -> ..temp.merged.dt;
            ##
            ## determine the final 'has.ge.1' in the union of original and edited
            ## NOTE: it is possible for the.motif.in.the.original.has.ge.1.site.overlapping.the.edit or the.motif.in.the.edited.has.ge.1.site.overlapping.the.edit has NA values. In this case we use the following rule: "[TRUE, NA] -> TRUE" and "[FALSE, NA] -> FALSE", because a 'NA' in original/edited set indicates that the motif has ==0 site overlapping this edit in that group, so the decision could be made based on the non-NA value solely. In pracice, it is equivalent to replace the NA's with FALSE first and then use the "|" operator directly, as "TRUE | FALSE == TRUE" and "FALSE | FALSE == FALSE"
            ..temp.merged.dt[is.na(the.motif.in.the.original.has.ge.1.site.overlapping.the.edit) == TRUE, the.motif.in.the.original.has.ge.1.site.overlapping.the.edit:=FALSE] -> ..temp.merged.dt
            ..temp.merged.dt[is.na(the.motif.in.the.edited.has.ge.1.site.overlapping.the.edit) == TRUE, the.motif.in.the.edited.has.ge.1.site.overlapping.the.edit:=FALSE] -> ..temp.merged.dt
            ..temp.merged.dt[
              , `:=`(
                    the.motif.in.the.union.of.original.and.edited.has.ge.1.site.overlapping.the.edit=(
                        the.motif.in.the.original.has.ge.1.site.overlapping.the.edit |
                        the.motif.in.the.edited.has.ge.1.site.overlapping.the.edit
                    )
                )
            ] ->  ..temp.merged.dt
            ##
            ## add transcript and edit info afterwards here so that we know which edit lost which miRNA binding sites
            data.table(
                TEMP.ROW.DT,
                ..temp.merged.dt
            ) -> TEMP.RESULT.DT
            ##
            TEMP.RESULT.DT
        } -> .; rbindlist(., use.names=TRUE) ->.;
        ##
        .
    } -> ..edited.motif.count.compared.with.original.dt
    ## here some transcripts might have zero detected motifs, because we scanned all whose 3'-UTR sequence can be put into the MAF
    ##
    ## confirm that there's no duplicated records for each (transcript.id, edit.POS, edit.rel.POS.wrt.3UTR, motif.name)
    ## > ..edited.motif.count.compared.with.original.dt[, .N, list(transcript.id, edit.POS, edit.rel.POS.wrt.3UTR, motif.name)][N>1]
    ## Empty data.table (0 rows and 5 cols): transcript.id,edit.POS,edit.rel.POS.wrt.3UTR,motif.name,N

    fwrite(..edited.motif.count.compared.with.original.dt, glue("{output.directory}/edited.motif.count.compared.with.original.dt.gz"))
    ..edited.motif.count.compared.with.original.dt
} -> edited.motif.count.compared.with.original.dt


{
    "./external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf" -> .;
    fread(glue('cat {.} | grep -v "#" |grep -P "\ttranscript\t" | sed -E \'s/.*gene_id "([^;]+)";.*transcript_id "([^;]+)";.*gene_name "([^;]+)".*/\\1\t\\2\t\\3/\' '), header=FALSE, col.names=c("gene.id", "transcript.id", "gene.name")) -> .;
    .
} -> gene.id.and.transcript.id.and.name.dt

merge(x=edited.motif.count.compared.with.original.dt, y=gene.id.and.transcript.id.and.name.dt, by="transcript.id", all.x=TRUE, all.y=FALSE) -> edited.motif.count.compared.with.original.with.gene.id.and.gene.name.dt

RE.valid.genes.only.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")

merge(x=edited.motif.count.compared.with.original.with.gene.id.and.gene.name.dt[the.motif.in.the.union.of.original.and.edited.has.ge.1.site.overlapping.the.edit==TRUE], y=RE.valid.genes.only.dt[is.na(AC)==FALSE, list(Gene_ID, POS, ID, SAMPLE, stage, is.normal, AF)], by.x=c("gene.id", "edit.POS"), by.y=c("Gene_ID", "POS"), all.x=TRUE, all.y=FALSE) -> edited.motif.count.compared.with.original.with.gene.id.and.gene.name.altered.only.RE.only.dt

fwrite(edited.motif.count.compared.with.original.with.gene.id.and.gene.name.altered.only.RE.only.dt, glue("{output.directory}/edited.motif.count.compared.with.original.with.gene.id.and.gene.name.altered.only.RE.only.dt.gz"))


{

    copy(edited.motif.count.compared.with.original.with.gene.id.and.gene.name.altered.only.RE.only.dt) -> .;
    .[, motif.site.count.change.upon.editing:=c("0->1"="gained", "1->0"="lost")[paste(sep="", motif.site.count.from.original, "->", motif.site.count.from.edited)] ] -> .;
    ## consider normal early stage samples only
    .[is.normal==TRUE][stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell", "8-cell")] -> .;
    ##= prettify stages
    temp.stage.dt <- data.table(read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties"))
    merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])] ->.;
    ##
    .[, motif.name.prettified:=c("PRE"="PRE", "ARE.9mer"="ARE\n(9-mer)", "Geissler2016"="hnRNPA1\nand A2/B1\nbinding site")[motif.name]] -> .;
    ##
    . -> ..to.plot.dt
    ##
    ggplot(..to.plot.dt, aes(x=paste(sep="", "<i>", gene.name, "</i><br/>", ID), y=AF, fill=motif.site.count.change.upon.editing)) -> .;
    . + geom_boxplot() ->.;
    . + coord_flip() -> .;
    . + facet_grid(motif.name.prettified~stage.description.ordered, scales="free_y", space="free_y") -> .;
    . + labs(x="", y="editing level", fill="") -> .;
    . + theme_pubr(base_size=10) ->.;
    . + theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_markdown()) -> .;
    . + scale_y_continuous(limits=c(0.1, 1), breaks=seq(0.1, 1, 0.3)) -> .;
    . + scale_fill_manual(values=c("#ED7D31", "#70AD47")) -> .;
    ggsave.A4(
        filename=glue("{output.directory}/RE.CCR4.NOT.motif.in.normal.early.stages.AF.boxplot.png"),
        plot=.,
        width.r=0.95, height.r=0.5
    )

    
    

    
}
