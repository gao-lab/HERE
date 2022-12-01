library("data.table")
library("magrittr")
source("./scripts/common/ggpubr.A4.R")

normal.edits.annotation.dt <- fread("./result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/snpEff.annotation.for.subset.observed.edits.dt.txt.gz")

normal.recoding.edits.dt <- {
    normal.edits.annotation.dt -> .;
    ##s select recoding edits
    .[grepl("missense|splice|stop_lost|start_lost|start_codon_gain", Annotation)] -> .;
}

fwrite(normal.recoding.edits.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/normal.recoding.edits.dt.csv.gz")

