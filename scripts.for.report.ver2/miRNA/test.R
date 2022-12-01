## Following the result of `run_internal_prepare_ts.R`

## 1. check gain
{
    ## pick changed combinations (transcript.id x edit.POS x edit.rel.POS.wrt.3UTR x miRNA.family.ID)
    edited.ts.human.compared.with.original.annotated.dt[miRNA.site.affected.type %in% c("gained")] -> .;
    ## merge with edited.ts.human.dt
    merge(x=edited.ts.human.dt, y=.,
          by=c("transcript.id", "edit.POS", "edit.rel.POS.wrt.3UTR", "miRNA.family.ID"),
          all.x=FALSE, all.y=TRUE) -> .;
}


## 2. check lost
