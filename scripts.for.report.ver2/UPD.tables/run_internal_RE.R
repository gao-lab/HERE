library("data.table")
library("magrittr")
source("./scripts/common/ggpubr.A4.R")

normal.RE.dt <- fread("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.dt.txt.gz")
UPD.recurrence.CJ.dt <- fread("./result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/subset.site.recurrence.comparison.CJ.dt.txt.gz")



normal.RE.edit.and.stage.only.dt <- {
    normal.RE.dt -> .;
    ##s select stages of interest only
    .[stage %in% c("zygote", "2-cell", "4-cell", "8-cell", "morula")] -> .;
    ##c keep edit x stage info only
    .[, list(CHROM, POS, stage, is.normal.RE=TRUE)] %>% unique -> .;
}

UPD.recurrence.CJ.pathological.only.dt <- {
    UPD.recurrence.CJ.dt -> .;
    ##s keep pathological samples
    .[grepl("androgenetic|parthenogenetic", group)] ->.;
    ##= extract stage and disease info
    .[, stage:=sub(pattern="@.*", replacement="", x=group)]
    .[, disease:=sub(pattern=".*@", replacement="", x=group)]
 }

UPD.disease.edit.and.normal.RE.comparison.dt <- merge(
    x=UPD.recurrence.CJ.pathological.only.dt, y=normal.RE.edit.and.stage.only.dt,
    by=c("CHROM", "POS", "stage"),
    all.x=TRUE, all.y=TRUE)[is.na(is.normal.RE)==TRUE, is.normal.RE:=FALSE]

fwrite(UPD.disease.edit.and.normal.RE.comparison.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/UPD.disease.edit.and.normal.RE.comparison.dt.csv.gz")

normal.RE.lost.in.UPD.dt <- UPD.disease.edit.and.normal.RE.comparison.dt[is.normal.RE==TRUE & recurrence.type != "recurrent"][, is.normal.RE:=NULL]

fwrite(normal.RE.lost.in.UPD.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/normal.RE.lost.in.UPD.dt.csv.gz")

