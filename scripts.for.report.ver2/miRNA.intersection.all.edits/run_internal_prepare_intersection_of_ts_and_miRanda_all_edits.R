library("data.table")
library("magrittr")
library("foreach")
library("doMC")
library("iterators")
library("glue")
library("stringi")
source("./scripts/common/ggpubr.A4.R")

output.directory <- "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA.intersection.all.edits/"
dir.create(output.directory, recursive=TRUE)

## 1. get original TargetScan results
{
    "result/A02_8__get_editing_effect_on_miRNA_binding_sites/step09__concatenate_TargetScan_results_across_all_chromosomes/32/gencode.3utr.all.chromosomes.concatenated.headless.TargetScan.output.gz" -> .;
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
        UTR.end=UTR_end,
        Site.type=Site_type)] -> .;
    ##= setkey
    setkey(., "transcript.id")
} -> original.TargetScan.human.dt
fwrite(original.TargetScan.human.dt, glue("{output.directory}/original.TargetScan.human.dt.gz"))

## 2. get edited TargetScan results (all edits)
{
    "result/A02_8__get_editing_effect_on_miRNA_binding_sites/step20__concatenate_all_edited_TargetScan_results_across_all_chromosomes/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/32/all.edited.gencode.3utr.all.chromosomes.concatenated.headless.TargetScan.output.gz" -> .;
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
        UTR.end=UTR_end,
        Site.type=Site_type)] ->.;
    ##= setkey
    setkey(., "transcript.id", "edit.POS")
} -> all.edited.TargetScan.human.dt
fwrite(all.edited.TargetScan.human.dt, glue("{output.directory}/all.edited.TargetScan.human.dt.gz"))


## "
## how each TargetScan site type is visualized in miRanda alignments 
## 1. site types with exact 7mer hits

## 1.1. 8mer-1a
##   [19, 26] is the hit region, with [19, 25] being the exact match to the seed region ([2, 7]) + position 8 of the microRNA, and [26] (NOTE: this is on the 3'-UTR, not on the microRNA) being the appended 'A' 
##      1             15  19     26
##   5' TAGCAGAGAGTCCTGAGCCACTGCCAACA 3'  ENST00000397026.7
##                        |||||||
##   3'                   GUGACGG     5'  miR-34-5p/449-5p
##                        8     2

## 1.2. 7mer-m8
##   [16, 22] is the hit region, with [16, 22] being the exact match to the seed region ([2, 7]) + position 8 of the microRNA, and [23] not being the appended 'A' (if [23] is 'A' then this hit will be 8mer-1a instead)
##      1              16    22
##   5' TGACAGGCACAGAGGGTGCCTTTTACCGCCGC 3' ENST00000456506.2
##                     |||||||
##   3'                CACGGAA           5' miR-124-3p.1
##                     8     2

## 2. site types with exact 6mer hits

## 2.1. 7mer-1a
##   [17, 23] is the hit region, with [17, 22] being the exact match to the seed region ([2, 7]) of the microRNA, and [23] being the appended 'A'
## ENST00000287652.8	miR-138-5p	9606	17	23	17	23	12957	7mer-1a	x	7mer-1a	9544 9606		0
##      1               17    23
##   5' TAAGCAGAGAAGGGCGACCAGCAGCCTGAC 3' ENST00000287652.8
##                      ||||||
##   3'                GUGGUCG         5' miR-138-5p
##                      7    2

## 2.2. 6mer  
##   [15, 20] is the hit region, with [15, 20] being the exact match to the seed region ([2, 7]) of the microRNA, and [21] not being the appended 'A' (if [21] is 'A' then this hit will be 7mer-1a instead)
## ENST00000425863.5	miR-24-3p	9606	15	20	15	20	2779	6mer	x	6mer	9544 9606		1
##      1             15   20
##   5' taataaaaatgcgttgagccta 3' ENST00000425863.5
##                    ||||||                
##   3'              GACUCGG   5' miR-24-3p
##                    7    2
 
## "

## function for 3 and 4

fread.miRanda.with.TargetScan.annotation <- function(filename){
    ## filename <- "result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step05__concatenate_miRanda_results_across_all_chromosomes/32/gencode.3utr.all.chromosomes.concatenated.headless.miRanda.output.gz"
    ##
    ## 1. load the raw miRanda table (with alignments)
    {
        print(glue("{date()} : reading table from {filename}"))
        filename -> .;
        fread(., header=FALSE, sep=",", select=c(1,2,3,4,5,10), col.names=c("miRNA.3to5.alignment", "alignment.details", "transcript.3UTR.alignment", "miRNA.name.and.miRNA.family.ID", "transcript.id.3UTR", "UTR.start.and.end"), strip.white=FALSE) -> .;
    } -> ..temp.raw.dt;
    ##
    ## 2. find 7mer sites (8mer-1a and 7mer-m8)
    {
        print(glue("{date()} : finding 7mer sites"))
        copy(..temp.raw.dt) -> .;
        ## 2.1. find likely 7mer hits (a contiguous range: [2, 8]=7mer in microRNA in the alignment)
        .[, stri_locate_first_regex(str=miRNA.3to5.alignment, pattern="[AGCU]{7}-*[AGCUagcu]{1}-*$")] -> ..TEMP.LIKELY.7MER.LOCATION.MATRIX;
        ## 2.2. find real 7mer hits (a likely 7mer hit, with all 7 bases in the 7mer matched to the transcript 3'-UTR)
        ## note that this is equivalent to searching the alignment detail for "^|||||||" IN THE PART OF the [2, 8] range only (the "^" must be present; otherwise it might search for [1, 7])
        .[, `:=`(likely.7mer.start=..TEMP.LIKELY.7MER.LOCATION.MATRIX[, "start"], likely.7mer.end=..TEMP.LIKELY.7MER.LOCATION.MATRIX[, "end"])] -> .;
        .[, likely.7mer.alignment.details := stri_sub(str=alignment.details, from=likely.7mer.start, to=likely.7mer.end)] -> .;
        .[stri_detect_regex(str=likely.7mer.alignment.details, pattern='^\\|{7}') == TRUE, contains.7mer:=TRUE][is.na(contains.7mer) == TRUE, contains.7mer:=FALSE] -> .;
        ## 2.3. classify 7mer hits into 8mer-1a and 7mer-m8
        .[, transcript.3UTR.likely.7mer.trailings:=stri_sub(str=transcript.3UTR.alignment, from=likely.7mer.start + 7 )] -> .;
        .[, transcript.3UTR.likely.7mer.first.nongap:=stri_match_first_regex(str=transcript.3UTR.likely.7mer.trailings, pattern="[AGCTagct]{1}")] -> .;
        .[(contains.7mer == TRUE) & (transcript.3UTR.likely.7mer.first.nongap %in% c("A", "a") == TRUE), site.type:="8mer-1a"] -> .; 
        .[(contains.7mer == TRUE) & (transcript.3UTR.likely.7mer.first.nongap %in% c("C", "G", "T", "c", "g", "t") == TRUE), site.type:="7mer-m8"] -> .;
        ##
        .
    } -> ..temp.with.7mer.sites.dt;
    ##
    ## 3. find 6mer sites (7mer-1a and 6mer)
    {
        print(glue("{date()} : finding 6mer sites"))
        copy(..temp.with.7mer.sites.dt) -> .;
        ## 3.1. find likely 6mer hits (a contiguous range: [2, 7]=6mer in microRNA in the alignment)
        .[, stri_locate_first_regex(str=miRNA.3to5.alignment, pattern="[AGCU]{6}-*[AGCUagcu]{1}-*$")] -> ..TEMP.LIKELY.6MER.LOCATION.MATRIX;
        ## 3.2. find real 6mer hits ( a likely 6mer hit, with all 6 bases in the 6mer matched to the transcript 3'-UTR
#### note that this is equivalent to searching the alignment detail for "^||||||" IN THE PART OF the [2, 7] range only (the "^" must be present; otherwise it might search for [1, 6])
        .[, `:=`(likely.6mer.start=..TEMP.LIKELY.6MER.LOCATION.MATRIX[, "start"], likely.6mer.end=..TEMP.LIKELY.6MER.LOCATION.MATRIX[, "end"])] -> .;
        .[, likely.6mer.alignment.details := stri_sub(str=alignment.details, from=likely.6mer.start, to=likely.6mer.end)] -> .;
        .[stri_detect_regex(str=likely.6mer.alignment.details, pattern='^\\|{6}') == TRUE, contains.6mer:=TRUE][is.na(contains.6mer) == TRUE, contains.6mer:=FALSE] -> .;
        ## 3.3. classify 6mer hits into 7mer-1a and 6mer
        .[, transcript.3UTR.likely.6mer.trailings:=stri_sub(str=transcript.3UTR.alignment, from=likely.6mer.start + 6 )] -> .;
        .[, transcript.3UTR.likely.6mer.first.nongap:=stri_match_first_regex(str=transcript.3UTR.likely.6mer.trailings, pattern="[AGCTagct]{1}")] -> .;
        .[(contains.7mer == FALSE) & (contains.6mer == TRUE) & (transcript.3UTR.likely.6mer.first.nongap %in% c("A", "a") == TRUE), site.type:="7mer-1a"] -> .; 
        .[(contains.7mer == FALSE) & (contains.6mer == TRUE) & (transcript.3UTR.likely.6mer.first.nongap %in% c("C", "G", "T", "c", "g", "t") == TRUE), site.type:="6mer"] -> .;
        ##
        .
    } -> ..temp.with.7mer.and.6mer.sites.dt
    ##
    ## 4. deduce the 3'-UTR site positions
    ## we need to delete gaps in the 3'-UTR when computing the positions
    ##
    ##    1033                 1052  reported UTR start and end
    ##    1033                 1054=1033+22-1 UTR start and end in the alignments (note that UTR is gapped)
    ##    
    ##    1    6    11  15     22
    ## 5' AGCTGGGT--TCTTCTACCTCA 3'
    ##    |:||: ::  :|| ||||||||
    ## 3' UUGAUAUGUUGGAUGAUGGAGU 5'
    ##    AGCTGGGT--TCTT            likely 7mer preceding segment
    {
        print(glue("{date()} : deducing 3'-UTR positions"))
        ## 4.1. get the UTR start
        copy(..temp.with.7mer.and.6mer.sites.dt) -> .;
        .[, UTR.start:=as.integer(sub(pattern="([0-9]+) ([0-9]+)", replacement="\\1", x=UTR.start.and.end))] -> .;
        ## 4.2. deduce the 3'-UTR site positions for 7mer sites
        .[contains.7mer == TRUE, likely.7mer.preceding.segment:=stri_sub(str=transcript.3UTR.alignment, from=1, to=likely.7mer.start - 1)] -> .;
        .[contains.7mer == TRUE, likely.7mer.preceding.segment.length:=stri_count_regex(str=likely.7mer.preceding.segment, pattern="[AGCTagct]{1}")] -> .;
        .[contains.7mer == TRUE, transcript.3UTR.site.start:= UTR.start + likely.7mer.preceding.segment.length] -> .;
        .[site.type == "8mer-1a", transcript.3UTR.site.end:=transcript.3UTR.site.start + 8 - 1] -> .;
        .[site.type == "7mer-m8", transcript.3UTR.site.end:=transcript.3UTR.site.start + 7 - 1] -> .;
        ## 4.3. deduce the 3'-UTR site positions for 6mer sites
        .[contains.7mer == FALSE & contains.6mer == TRUE, likely.6mer.preceding.segment:=stri_sub(str=transcript.3UTR.alignment, from=1, to=likely.6mer.start - 1)] -> .;
        .[contains.7mer == FALSE & contains.6mer == TRUE, likely.6mer.preceding.segment.length:=stri_count_regex(str=likely.6mer.preceding.segment, pattern="[AGCTagct]{1}")] -> .;
        .[contains.7mer == FALSE & contains.6mer == TRUE, transcript.3UTR.site.start:= UTR.start + likely.6mer.preceding.segment.length] -> .;
        .[site.type == "7mer-1a", transcript.3UTR.site.end:=transcript.3UTR.site.start + 7 - 1] -> .;
        .[site.type == "6mer", transcript.3UTR.site.end:=transcript.3UTR.site.start + 6 - 1] -> .;
        ##
        .
    } -> ..temp.with.7mer.and.6mer.sites.with.3UTR.positions.dt;
    print(glue("{date()} : finished"))
    ## 
    return(..temp.with.7mer.and.6mer.sites.with.3UTR.positions.dt)
}

## 3. get original miRanda results
{
    ##
    "result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step05__concatenate_miRanda_results_across_all_chromosomes/32/gencode.3utr.all.chromosomes.concatenated.headless.miRanda.output.gz" -> .;
    fread.miRanda.with.TargetScan.annotation(.) -> .;
    ##
    ## formalize columns
    .[, list(
        transcript.id=sub(pattern="__3UTR$", replacement="", x=transcript.id.3UTR),
        miRNA.name=sub(pattern="__.*", replacement="", x=miRNA.name.and.miRNA.family.ID),
        miRNA.family.ID=sub(pattern=".*__", replacement="", x=miRNA.name.and.miRNA.family.ID),
        UTR.start=UTR.start,
        UTR.end=as.integer(sub(pattern="([0-9]+) ([0-9]+)", replacement="\\2", x=UTR.start.and.end)),
        transcript.3UTR.site.start=transcript.3UTR.site.start,
        transcript.3UTR.site.end=transcript.3UTR.site.end,
        site.type=site.type
    )] -> .;
    setkey(., "transcript.id")
    ##
    .
} -> original.miRanda.human.dt
fwrite(original.miRanda.human.dt, glue("{output.directory}/original.miRanda.human.dt.gz"))

## 4. get edited miRanda results
{
    ##
    "result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step16__concatenate_all_edited_miRanda_results_across_all_chromosomes/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/32/all.edited.gencode.3utr.all.chromosomes.concatenated.headless.miRanda.output.gz" -> .;
    fread.miRanda.with.TargetScan.annotation(.) -> .;
    ##
    ##= extract columns and formalize column names
    .[, list(
        transcript.id=sub(pattern="^([^@]+)@.*", replacement="\\1", x=transcript.id.3UTR),
        edit.POS=sub(pattern="^[^@]+@([^@]+)@.*", replacement="\\1", x=transcript.id.3UTR) %>% as.integer,
        edit.rel.POS.wrt.3UTR=sub(pattern="^[^@]+@[^@]+@([^@]+)__3UTR$", replacement="\\1", x=transcript.id.3UTR) %>% as.integer,
        miRNA.name=sub(pattern="__.*", replacement="", x=miRNA.name.and.miRNA.family.ID),
        miRNA.family.ID=sub(pattern=".*__", replacement="", x=miRNA.name.and.miRNA.family.ID),
        UTR.start=UTR.start,
        UTR.end=as.integer(sub(pattern="([0-9]+) ([0-9]+)", replacement="\\2", x=UTR.start.and.end)),
        transcript.3UTR.site.start=transcript.3UTR.site.start,
        transcript.3UTR.site.end=transcript.3UTR.site.end,
        site.type=site.type
    )] -> .;
    ##= setkey
    setkey(., "transcript.id", "edit.POS")
    ##
    .
} -> all.edited.miRanda.human.dt
fwrite(all.edited.miRanda.human.dt, glue("{output.directory}/all.edited.miRanda.human.dt.gz"))


## 5. merge prediction, original TargetScan + original miRanda
## only merge those with same site type and same 3UTR region
## see case.validation.org for case validation
{
    ## [1] merge the two tables
    {
        merge(
            x=original.TargetScan.human.dt[, list(transcript.id, miRNA.family.ID, transcript.3UTR.site.start=UTR.start, transcript.3UTR.site.end=UTR.end, site.type=Site.type, predicted.by.TargetScan=TRUE, prediction.index.TargetScan=.I)],
            y=original.miRanda.human.dt[, list(transcript.id, miRNA.family.ID, transcript.3UTR.site.start, transcript.3UTR.site.end, site.type, predicted.by.miRanda=TRUE, prediction.index.miRanda=.I, miRNA.name.from.miRanda.input=miRNA.name, miRanda.UTR.start=UTR.start, miRanda.UTR.end=UTR.end)],
            by=c("transcript.id", "miRNA.family.ID", "transcript.3UTR.site.start", "transcript.3UTR.site.end", "site.type"),
            all=TRUE, allow.cartesian=FALSE) -> .;
        ## all combinations have no duplicates
        ## .[, .N, list(transcript.id, miRNA.family.ID, transcript.3UTR.site.start, transcript.3UTR.site.end, site.type, miRNA.name.from.miRanda.input, miRanda.UTR.start, miRanda.UTR.end)][N>1]
        ## Empty data.table (0 rows and 9 cols): transcript.id,miRNA.family.ID,transcript.3UTR.site.start,transcript.3UTR.site.end,site.type,miRNA.name.from.miRanda.input...
        .
    } -> ..original.TargetScan.and.miRanda.match.dt
    ## [1].backup
    fwrite(..original.TargetScan.and.miRanda.match.dt, glue("{output.directory}/original.TargetScan.and.miRanda.match.dt.gz"))
    ##
    ## [2] take the intersection
    ## take the intersection
    ..original.intersection.of.TargetScan.and.miRanda.dt <- ..original.TargetScan.and.miRanda.match.dt[(predicted.by.TargetScan == TRUE) & (predicted.by.miRanda == TRUE)]
    ## [2].backup
    fwrite(..original.intersection.of.TargetScan.and.miRanda.dt, glue("{output.directory}/original.intersecion.of.TargetScan.and.miRanda.match.dt.gz"))
    ## [3] collapse to miRNA family level (this is important for evaluating the change of MBS per miRNA family, because while TargetScan can only predict for each miRNA family rather than for each individual mature miRNA, miRanda can predict for each individual mature miRNA; as a result, a single prediction in TargetScan might match multiple mature miRNA prediction in miRanda, whose site regions are all identical to that in TargetScan. A fair match should be done at the miRNA family level)
    {
        copy(..original.intersection.of.TargetScan.and.miRanda.dt) -> .;
        .[,
          list(
              prediction.indices.miRanda=paste(collapse=";", prediction.index.miRanda),
              all.miRNA.names.from.miRanda.input=paste(collapse=";", miRNA.name.from.miRanda.input),
              all.miRanda.UTR.starts.and.ends=paste(collapse=";", paste(sep="", miRanda.UTR.start, "-", miRanda.UTR.end))
          ),
          list(
              transcript.id, 
              miRNA.family.ID,
              transcript.3UTR.site.start, transcript.3UTR.site.end, site.type,
              predicted.by.TargetScan, prediction.index.TargetScan,
              predicted.by.miRanda
        )] -> .;
    } -> ..original.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt
    ## [3].backup
    fwrite(..original.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt, glue("{output.directory}/original.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt.gz"))
    ..original.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt
} -> original.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt


## 6. merge prediction, edited TargetScan + edited miRanda
## only merge those with same site type and same 3UTR region
{
    ## [1] merge two tables
    {
        merge(
            x=all.edited.TargetScan.human.dt[, list(transcript.id, edit.POS, edit.rel.POS.wrt.3UTR, miRNA.family.ID, transcript.3UTR.site.start=UTR.start, transcript.3UTR.site.end=UTR.end, site.type=Site.type, predicted.by.TargetScan=TRUE, prediction.index.TargetScan=.I)],
            y=all.edited.miRanda.human.dt[, list(transcript.id, edit.POS, edit.rel.POS.wrt.3UTR, miRNA.family.ID, transcript.3UTR.site.start, transcript.3UTR.site.end, site.type, predicted.by.miRanda=TRUE, prediction.index.miRanda=.I, miRNA.name.from.miRanda.input=miRNA.name, miRanda.UTR.start=UTR.start, miRanda.UTR.end=UTR.end)],
            by=c("transcript.id", "edit.POS", "edit.rel.POS.wrt.3UTR", "miRNA.family.ID", "transcript.3UTR.site.start", "transcript.3UTR.site.end", "site.type"),
            all=TRUE, allow.cartesian=FALSE) -> .;
        ## all combinations have no duplicates
        ## .[, .N, list(transcript.id, edit.POS, edit.rel.POS.wrt.3UTR,  miRNA.family.ID, transcript.3UTR.site.start, transcript.3UTR.site.end, site.type, miRNA.name.from.miRanda.input, miRanda.UTR.start, miRanda.UTR.end)][N>1]
        ## Empty data.table (0 rows and 11 cols): transcript.id,edit.POS,edit.rel.POS.wrt.3UTR,miRNA.family.ID,transcript.3UTR.site.start,transcript.3UTR.site.end...
    } -> ..all.edited.TargetScan.and.miRanda.match.dt
    ## [1].backup
    fwrite(..all.edited.TargetScan.and.miRanda.match.dt, glue("{output.directory}/all.edited.TargetScan.and.miRanda.match.dt.gz"))
    ##
    ## [2] annotate whether each edit is on each site
    {
        ..all.edited.TargetScan.and.miRanda.match.dt[, edit.is.within.site:=(
            (edit.rel.POS.wrt.3UTR >= transcript.3UTR.site.start) &
            (edit.rel.POS.wrt.3UTR <= transcript.3UTR.site.end)
        )] -> .;
    } -> ..all.edited.TargetScan.and.miRanda.match.with.edit.site.residence.dt
    ## [3] take the intersection
    {
        ..all.edited.TargetScan.and.miRanda.match.dt[(predicted.by.TargetScan == TRUE) & (predicted.by.miRanda == TRUE)] -> .;
    } -> ..all.edited.intersection.of.TargetScan.and.miRanda.dt
    ## [3].backup
    fwrite(..all.edited.intersection.of.TargetScan.and.miRanda.dt, glue("{output.directory}/all.edited.intersecion.of.TargetScan.and.miRanda.match.dt.gz"))
    ## [4] collapse to the miRNA family level
    ## the edited version differs from the original version in the columns (it used 3 additional columns, edit.POS and edit.rel.POS.wrt.3UTR, and edit.is.within.site)
    {
        copy(..all.edited.intersection.of.TargetScan.and.miRanda.dt) -> .;
        .[,
          list(
              prediction.indices.miRanda=paste(collapse=";", prediction.index.miRanda),
              all.miRNA.names.from.miRanda.input=paste(collapse=";", miRNA.name.from.miRanda.input),
              all.miRanda.UTR.starts.and.ends=paste(collapse=";", paste(sep="", miRanda.UTR.start, "-", miRanda.UTR.end)),
              the.site.in.the.edited.has.ge.1.site.overlapping.the.edit=any(edit.is.within.site) ## mark whether this miRNA family has any sites overlapping the edit (for piechart)
          ),
          list(
              transcript.id, edit.POS, edit.rel.POS.wrt.3UTR,
              miRNA.family.ID,
              transcript.3UTR.site.start, transcript.3UTR.site.end, site.type,
              predicted.by.TargetScan, prediction.index.TargetScan,
              predicted.by.miRanda
        )] -> .;
    } -> ..all.edited.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt
    ## [4].backup
    fwrite(..all.edited.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt, glue("{output.directory}/all.edited.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt.gz"))
    ## 
    ..all.edited.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt
} -> all.edited.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt



## 8. merge the two to get a list of gain/lost miRNA binding sites per transcript x edit POS
{
    
    ## [1] get all (transcript.id, edit.POS) pairs
    {
        ##
        "./result/A02_8__get_editing_effect_on_miRNA_binding_sites/step16__compute_all_edit_relative_position_on_3UTR/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/32/gencode.3utr.and.all.edit.CJ.dt.csv.gz" -> .;
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
    } -> ..all.transcript.and.all.edit.pairs.to.check.dt
    
    ## [2] base table (transcript x edit)    
    {
        ##
        registerDoMC(cores=7)
        foreach(TEMP.ROW.DT=iter(..all.transcript.and.all.edit.pairs.to.check.dt, by="row")) %dopar% {
            TEMP.INDEX <- TEMP.ROW.DT[1, index]
            if (TEMP.INDEX %% 1000 == 0){
                print(glue("{date()} : processing index {TEMP.INDEX}"))
            }
            ##
            ## get original hit counts (in two steps rather than one step, because we need to annotate for each miRNA family the.site.in.the.original.has.ge.1.site.overlapping.the.edit)
            original.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt[
                transcript.id == TEMP.ROW.DT[1, transcript.id]
            ][,
              the.site.in.the.original.has.ge.1.site.overlapping.the.edit:=(
                  (TEMP.ROW.DT[1, edit.rel.POS.wrt.3UTR] >= transcript.3UTR.site.start) &
                  (TEMP.ROW.DT[1, edit.rel.POS.wrt.3UTR] <= transcript.3UTR.site.end)
              )] -> ..temp.original.hit.dt
            ##
            ..temp.original.hit.dt[
              , list(
                    miRNA.family.site.count.from.original=.N,
                    from.original=TRUE,
                    the.miRNA.family.in.the.original.has.ge.1.site.overlapping.the.edit=any(the.site.in.the.original.has.ge.1.site.overlapping.the.edit)),
                list(miRNA.family.ID)
            ] -> ..temp.original.hit.count.dt
            ##
            ## get edited hit counts
            all.edited.intersection.of.TargetScan.and.miRanda.at.miRNA.family.level.dt[
                (transcript.id == TEMP.ROW.DT[1, transcript.id]) & (edit.POS == TEMP.ROW.DT[1, edit.POS]),
                list(
                    miRNA.family.site.count.from.edited=.N,
                    from.edited=TRUE,
                    the.miRNA.family.in.the.edited.has.ge.1.site.overlapping.the.edit=any(the.site.in.the.edited.has.ge.1.site.overlapping.the.edit)),
                list(miRNA.family.ID)
            ][is.na(miRNA.family.ID) == FALSE] -> ..temp.edited.hit.count.dt ## automatically create zero-row table if not in the edited MBS set
            ##
            ## merge all columns but transcript and edit info first
            merge(
                x=..temp.original.hit.count.dt, y=..temp.edited.hit.count.dt,
                by=c("miRNA.family.ID"),
                all.x=TRUE, all.y=TRUE
            ) -> ..temp.merged.dt
            ##
            ## fill the zeroes and FALSE's
            ..temp.merged.dt[is.na(from.original)==TRUE, `:=`(miRNA.family.site.count.from.original=0, from.original=FALSE)] -> ..temp.merged.dt;
            ..temp.merged.dt[is.na(from.edited)==TRUE, `:=`(miRNA.family.site.count.from.edited=0, from.edited=FALSE)] -> ..temp.merged.dt;
            ## determine the final 'has.ge.1' in the union of original and edited
            ## NOTE: it is possible for the.miRNA.family.in.the.original.has.ge.1.site.overlapping.the.edit or the.miRNA.family.in.the.edited.has.ge.1.site.overlapping.the.edit has NA values. In this case we use the following rule: "[TRUE, NA] -> TRUE" and "[FALSE, NA] -> FALSE", because a 'NA' in original/edited set indicates that the miRNA family has ==0 site overlapping this edit in that group, so the decision could be made based on the non-NA value solely. In pracice, it is equivalent to replace the NA's with FALSE first and then use the "|" operator directly, as "TRUE | FALSE == TRUE" and "FALSE | FALSE == FALSE"
            ..temp.merged.dt[is.na(the.miRNA.family.in.the.original.has.ge.1.site.overlapping.the.edit) == TRUE, the.miRNA.family.in.the.original.has.ge.1.site.overlapping.the.edit:=FALSE] -> ..temp.merged.dt
            ..temp.merged.dt[is.na(the.miRNA.family.in.the.edited.has.ge.1.site.overlapping.the.edit) == TRUE, the.miRNA.family.in.the.edited.has.ge.1.site.overlapping.the.edit:=FALSE] -> ..temp.merged.dt
            ..temp.merged.dt[
              , `:=`(
                    the.miRNA.family.in.the.union.of.original.and.edited.has.ge.1.site.overlapping.the.edit=(
                        the.miRNA.family.in.the.original.has.ge.1.site.overlapping.the.edit |
                        the.miRNA.family.in.the.edited.has.ge.1.site.overlapping.the.edit
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
    } -> ..all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.dt
    fwrite(..all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.dt, glue("{output.directory}/all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.dt.gz"))
    
    ## confirm that there's no transcript with zero detected families (those with zero families have been discarded in the loop above)
    ## ..edited.intersection.of.TargetScan.and.miRanda.dt[transcript.id %in% edited.TargetScan.and.miRanda.compared.with.original.dt[, list(count.of.edit.detected.families=sum(from.edited)), transcript.id][count.of.edit.detected.families==0, transcript.id]]
    ## Empty data.table (0 rows and 14 cols): transcript.id,edit.POS,edit.rel.POS.wrt.3UTR,miRNA.family.ID,transcript.3UTR.site.start,transcript.3UTR.site.end...
    ##
    ## confirm that there's no duplicated records for each (transcript.id, edit.POS, edit.rel.POS.wrt.3UTR, miRNA.family.ID)
    ## ..edited.intersection.of.TargetScan.and.miRanda.compared.with.original.dt[, .N, list(transcript.id, edit.POS, edit.rel.POS.wrt.3UTR, miRNA.family.ID)][N>1]
    ## Empty data.table (0 rows and 5 cols): transcript.id,edit.POS,edit.rel.POS.wrt.3UTR,miRNA.family.ID,N
    ..all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.dt
} -> all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.dt


## see case.validation.has.ge.1.org to check whether the the.miRNA.family.has.ge.1.site.overlapping.the.edit is correct by case validation (one positive and one negative case)

## 9. get annotated table (transcript x edit, with annotation)
{
    
    ## [1] annotate affected type
    {
        copy(all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.dt) -> .;
        ## remove is.na(miRNA.family.ID) (i.e., no MBS before, AND no MBS after editing)
        .[is.na(miRNA.family.ID) == FALSE] -> .;
        ## add miRNA site affected type
        .[miRNA.family.site.count.from.original == miRNA.family.site.count.from.edited, miRNA.site.affected.type:="unchanged"]
        .[miRNA.family.site.count.from.original < miRNA.family.site.count.from.edited, miRNA.site.affected.type:="gained"]
        .[miRNA.family.site.count.from.original > miRNA.family.site.count.from.edited, miRNA.site.affected.type:="lost"]
        ## add miRNA family ID affected type
        .[miRNA.family.site.count.from.original > 0 & miRNA.family.site.count.from.edited > 0, miRNA.family.ID.affected.type:="unchanged"]
        .[miRNA.family.site.count.from.original == 0 & miRNA.family.site.count.from.edited > 0, miRNA.family.ID.affected.type:="gained"]
        .[miRNA.family.site.count.from.original > 0 & miRNA.family.site.count.from.edited == 0, miRNA.family.ID.affected.type:="lost"]
        ##
        .
    } -> ..all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.dt
    ## [1].backup
    fwrite(..all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.dt, glue("{output.directory}/all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.dt.gz"))

    ## [2] summarize table (transcript x edit, with net effect of each edit on each transcript)
    {
        copy(..all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.dt) -> .;
        .[,
          list(
              difference.of.count.of.miRNA.sites=sum(miRNA.family.site.count.from.edited - miRNA.family.site.count.from.original),
              difference.of.count.of.miRNA.families=sum(miRNA.family.ID.affected.type=='gained') - sum(miRNA.family.ID.affected.type=='lost')
          ),
          list(transcript.id, edit.POS)
          ] -> .; 
    } -> ..all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.summary.dt
    ## [2].backup
    fwrite(..all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.summary.dt, glue("{output.directory}/all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.summary.dt.gz"))

    ## [3-pre] get gene-transcript table
    {
        "external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf" -> .;
        fread(cmd=glue("cat {.} | grep -v '^#' | grep -P '\ttranscript\t' | cut -f 9 | sed -E -e 's@.*gene_id \"([^\"]+)\"; transcript_id \"([^\"]+)\";.*@\\1\t\\2@' "), header=FALSE, sep="\t", col.names=c("gene.id", "transcript.id")) -> .;
        .
    } -> ..gene.transcript.mapping.dt

    ## [3] summarize table at gene level (gene x edit, with net effect of each edit on each gene)
    {
        copy(..all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.summary.dt) -> .;
        ##m merge with gene info
        merge(x=., y=..gene.transcript.mapping.dt,
              by.x="transcript.id", by.y="transcript.id",
              all.x=TRUE, all.y=FALSE) -> .;
        ##c collapse at the level of gene x edit
        .[,
          list(
              gains.miRNA.sites=all(difference.of.count.of.miRNA.sites >= 0) & any(difference.of.count.of.miRNA.sites > 0),
              gains.miRNA.families=all(difference.of.count.of.miRNA.families >= 0) & any(difference.of.count.of.miRNA.families > 0),
              loses.miRNA.sites=all(difference.of.count.of.miRNA.sites <= 0) & any(difference.of.count.of.miRNA.sites < 0),
              loses.miRNA.families=all(difference.of.count.of.miRNA.families <= 0) & any(difference.of.count.of.miRNA.families < 0)
          ),
          list(gene.id, edit.POS)
          ] -> .;
        ##
        .
    } -> ..all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.summary.gene.and.edit.level.dt
    ## [3].backup
    fwrite(..all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.summary.gene.and.edit.level.dt, glue("{output.directory}/all.edited.intersection.of.TargetScan.and.miRanda.compared.with.original.annotated.summary.gene.and.edit.level.dt.gz"))
    
}
