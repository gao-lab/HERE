## Following the result of `run_internal_prepare_ts.R` (till #### 3.2. annotated table (transcript x edit, with annotation); we need the set of 'gain')
## all coordinates (edit position, edit pos rel to xx, and MSA/UTR start and ends) are all one-based
## UTR counts from stop codon, with the first base of stop codon being coordinate 1

gained.dt <- edited.ts.human.compared.with.original.annotated.dt[miRNA.site.affected.type=='gained']

original.ts.full.human.dt <- {
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
        UTR.end=UTR_end,
        Group_type, Species_in_this_group, Species_in_this_group_with_this_site_type, ORF_overlap)] -> .;
    ##= setkey
    setkey(., "transcript.id")
}



edited.ts.full.human.dt <- {
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
        UTR.end=UTR_end,
        Group_type, Species_in_this_group, Species_in_this_group_with_this_site_type, ORF_overlap)] ->.;
    ##= setkey
    setkey(., "transcript.id", "edit.POS")
}


aa <- {
    copy(edited.ts.full.human.dt) -> .;
    ## get match in original
    .[paste(sep="", transcript.id, "@", miRNA.family.ID, "@", MSA.start, "@", MSA.end) %in% original.ts.human.dt[, paste(sep="", transcript.id, "@", miRNA.family.ID, "@", MSA.start, "@", MSA.end)], detected.in.original:=TRUE]
    .[is.na(detected.in.original) == TRUE, detected.in.original:=FALSE]
}

{
    aa[, `:=`(edit.rel.POS.wrt.MSA.start=edit.rel.POS.wrt.3UTR - MSA.start + 1)]
    aa[paste(sep="", transcript.id, "@", edit.POS, "@", miRNA.family.ID) %in% gained.dt[, paste(sep="", transcript.id, "@", edit.POS, "@", miRNA.family.ID)], `:=`(is.a.gain=TRUE)]
    aa[is.na(is.a.gain)==TRUE, is.a.gain:=FALSE]
    ## for some positive cases:
    aa[detected.in.original==FALSE][edit.rel.POS.wrt.MSA.start==2][MSA.start==29]
    ## for problematic` cases (but not a gain):
    aa[detected.in.original==FALSE][edit.rel.POS.wrt.MSA.start==0][MSA.start==28]
    ## for problematic` cases (and is a gain):
    aa[detected.in.original==FALSE][is.a.gain==TRUE][, .N, list(edit.rel.POS.wrt.MSA.start.is.ge.1=edit.rel.POS.wrt.MSA.start>=1, edit.rel.POS.wrt.MSA.start.is.le.UTR.end.rel.POS.wrt.UTR.start=edit.rel.POS.wrt.MSA.start <= (UTR.end-UTR.start+1))]
'
   edit.rel.POS.wrt.MSA.start.is.ge.1
1:                               TRUE
   edit.rel.POS.wrt.MSA.start.is.le.UTR.end.rel.POS.wrt.UTR.start    N
1:                                                           TRUE 6374

'
    ##TRUE x TRUE only
    ##validation passed: all new MBSs (with a strictly positive net increase in miRNA family; see the `problematic` case below) can cover the edit site

}


'
## positive cases (positive strand):
chr19:ENST00000451849.1 (positive strand, first base of stop codon: 57756877)
miR-130-3p/301-3p/454-3p (6mer)
miR-148-3p/152-3p (7mer-m8)
   57756905 (absolute position)
   29    35 (UTR start and end rel. to this transcript)
    *       (edit position)
  +TACACTG  (DNA, original strand, 5'-to-3')
   UACACUG  (unedited RNA, 5'-to-3')
   UGCACUG  (edited RNA, 5'-to-3')
~~~the following should not be directly aligned to the above, because the transcript is reversed~~~
   GUCACGU (edited RNA, 3'-to-5'(reverse))
   CAGUGCA (edited RNA, complementary strand, 3'-to-5'(reverse))
   _AGUGCAA (miRNA seed region for miR-130-3p/301-3p/454-3p (6mer), original strand, 5'-to-3')
   CAGUGCA (miRNA seed region for miR-130-3p/301-3p/454-3p (6mer), original strand, 5'-to-3')

## positive cases (negative strand):
chr22:ENST00000399568.5 (negative strand)
miR-150-5p
        19441622
   30   25 (UTR start and end rel. to this transcript)
       *  (edit position)
  +CTCCTA (DNA, original strand, 5'-to-3')
  -GAGGAT (DNA, complementary strand, 3'-to-5' (reverse))
   GAGGAU (unedited RNA, 3'-to-5'(reverse))
   GAGGGU (edited RNA, 3'-to-5'(reverse))
   CUCCCA (edited RNA, complementary strand, 3'-to-5'(reverse))
   CUCCCAA (miRNA seed region, original strand, 5'-to-3')


## NOT a problematic cases (rel.POS.wrt.MSA.start == 0, but the number of the same microRNA family site does not change upon editing -- a binding site in UTR 27-33 exists before editing, and while this site is lost after editing, a new binding site for the same micrRNA family ID emergs in UTR 28-33 after editing):
chr2:ENST00000481670.2 (negative strand, first base of stop codon: 58213628)
miR-338-3p
        58213627
   33   28  (UTR start and end rel. to this transcript)
         *  (edit position (27 here))
  +CCAGCAT (DNA, original strand, 5'-to-3')
  -GGTCGTA (DNA, complementary strand, 3'-to-5' (reverse))
   GGUCGUA (unedited RNA, 3'-to-5'(reverse))
   GGUCGUG (edited RNA, 3'-to-5'(reverse))
   CCAGCAC (edited RNA, complementary strand, 3'-to-5'(reverse))
   CCAGCAU  (miRNA seed region, original strand, 5'-to-3')

'
