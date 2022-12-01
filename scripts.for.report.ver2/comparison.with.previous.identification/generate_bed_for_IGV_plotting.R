library("data.table")


{
    
    "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/qiu2016.bed.merged.with.ours.dt.gz" -> .;
    fread(.) -> .;
    .[is.na(is.qiu2016) == TRUE, is.qiu2016:=FALSE] -> .;
    .[is.na(is.ours) == TRUE, is.ours:=FALSE] -> .;
    unique(.[, list(Chr.hg38, pos, is.qiu2016, is.ours)]) -> .;
    ## examine the TTF1 recoding site region
    .[Chr.hg38=='chr9'][pos >= 132375956-40 & pos <= 132375956+40] -> .;
    . -> info.dt
    
    ## generate bed files
    info.dt[is.qiu2016 == TRUE] -> .;
    .[, list(chrom=Chr.hg38, chromStart=pos-1, chromEnd=pos-1+1, name=".", score=500, strand=".")] -> .;
    . -> qiu2016.bed.dt
    fwrite(qiu2016.bed.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/case.TTF1.qiu2016.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

    ## generate bed files
    info.dt[is.ours == TRUE] -> .;
    .[, list(chrom=Chr.hg38, chromStart=pos-1, chromEnd=pos-1+1, name=".", score=500, strand=".")] -> .;
    . -> ours.bed.dt
    fwrite(ours.bed.dt, "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/total.samples/case.TTF1.ours.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

    ## plot region: chr9:132,375,916-132,375,996
}
