library("data.table")

dir.create("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA.case.study", recursive=TRUE)

input.reference.GTF.filename <- "external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf"
gene.id.and.symbol.mapping.dt <- fread(cmd=paste(sep="", "cat ", input.reference.GTF.filename, " | grep -v '^#' | grep -P '\tgene\t' | cut -f 9 | sed -E -e 's@.*gene_id \"([^\"]+)\";.*gene_name \"([^\"]+)\";.*@\\1\t\\2@' "), header=FALSE, sep="\t", col.names=c("gene.id", "gene.symbol"))

MTC.info.dt <- fread("result/S42_1__annotate_embryonic_genes/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/201221-fifth-phenotype-collection/combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt.txt.gz")[, list(gene.id, cluster)]


all.cases.dt <- {
    fread("./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/edited.ts.human.compared.with.original.annotated.summary.gene.and.edit.level.dt.csv.gz") ->.;
    merge(x=., y=gene.id.and.symbol.mapping.dt,
          by="gene.id",
          all.x=TRUE, all.y=FALSE) ->.;
    merge(x=., y=MTC.info.dt,
          by="gene.id",
          all.x=TRUE, all.y=FALSE)
}

{
    copy(all.cases.dt) ->.;
    ##s select SUV39H2
    .[gene.symbol=="SUV39H2"] ->.;
    ##= annotate type
    .[gains.miRNA.sites==TRUE & loses.miRNA.sites==FALSE, type:="MBS-gaining"]
    .[gains.miRNA.sites==FALSE & loses.miRNA.sites==TRUE, type:="MBS-losing"]
    .[gains.miRNA.sites==FALSE & loses.miRNA.sites==FALSE, type:="MBS-neutral"]
    .[, type.color:=c("MBS-gaining"="237,125,49", "MBS-losing"="112,173,71", "MBS-neutral"="68,114,196")[type]]
    ##s discard records without type
    .[is.na(type)==FALSE] ->.;
    ## generate bed-formatted record
    .[, list(CHROM="chr10", start=edit.POS-1, end=edit.POS-1+1, name=type, score=1000, strand=".", thickStart=edit.POS-1, thickEnd=edit.POS-1+1, itemRgb=type.color)] ->.;
    ## write to disk
    fwrite(., "./report.ver2/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/RE.miRNA/case.study.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}
