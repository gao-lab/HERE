library("data.table")
library("magrittr")

human.pc.gene.info.dt <- fread(cmd=paste(sep="", ' grep -v "#" ', snakemake@input[["reference_GTF_filename"]], ' | grep -P "\tgene\t" | grep -P \'gene_type "protein_coding"\' | cut -f 9 | sed -E -e \'s@gene_id "([^;]+)";.*gene_name "([^;]+)";.*@\\1\\t\\2@\' '), header=FALSE, col.names=c("gene.id", "gene.name"))

if (FALSE) {
    human.pc.gene.info.dt <- fread(cmd=paste(sep="", ' grep -v "#" ', 'external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf', ' | grep -P "\tgene\t" | grep -P \'gene_type "protein_coding"\' | cut -f 9 | sed -E -e \'s@gene_id "([^;]+)";.*gene_name "([^;]+)";.*@\\1\\t\\2@\' '), header=FALSE, col.names=c("gene.id", "gene.name") )
}


## combined.gexpr.FPKM.transposed.pc.only.dt <- snakemake@input[["combined_gexpr_FPKM_matrix_filename"]] %>%
##     {
##         temp.dt <- fread(., header=FALSE, skip=1, drop=1) %>%
##             as.matrix %>% t %>% as.data.table
##         setnames(temp.dt, fread(., header=FALSE, skip=1, select=1, col.names="gene.id")[, gene.id])
##         temp.pc.dt <- temp.dt[gene.id %in% human.pc.gene.info.dt[, gene.id]]
##         temp.pc.dt[, SAMPLE:=scan(file=., what=character(), nlines=1)]
##         temp.pc.dt
##     }

## if (FALSE) {

##     combined.gexpr.FPKM.transposed.pc.only.dt <- "result/BS06_1__get_expression_level/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/combined.gexpr.FPKM.matrix.txt" %>%
##     {
##         temp.dt <- fread(., header=FALSE, skip=1, drop=1) %>%
##             as.matrix %>% t %>% as.data.table
##         setnames(temp.dt, fread(., header=FALSE, skip=1, select=1, col.names="gene.id")[, gene.id])
##         temp.pc.dt <- temp.dt[, human.pc.gene.info.dt[, gene.id], with=FALSE]
##         temp.pc.dt[, SAMPLE:=scan(file=., what=character(), nlines=1)]
##         temp.pc.dt
##     }

## }

combined.gexpr.FPKM.pc.only.melt.dt <- snakemake@input[["combined_gexpr_FPKM_matrix_filename"]] %>%
    {
        temp.dt <- fread(., header=FALSE, skip=1)
        setnames(temp.dt, c("gene.id", scan(file=., what=character(), nlines=1)))
        temp.pc.dt <- temp.dt[gene.id %in% human.pc.gene.info.dt[, gene.id]]
        temp.pc.dt
    } %>%
    {melt(., id.vars="gene.id", variable.name="SAMPLE", value.name="FPKM")}


if (FALSE) {
    combined.gexpr.FPKM.pc.only.melt.dt <- "./result/BS06_1__get_expression_level/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/combined.gexpr.FPKM.matrix.txt" %>%
        {
            temp.dt <- fread(., header=FALSE, skip=1)
            setnames(temp.dt, c("gene.id", scan(file=., what=character(), nlines=1)))
            temp.pc.dt <- temp.dt[gene.id %in% human.pc.gene.info.dt[, gene.id]]
            temp.pc.dt
        } %>%
        {melt(., id.vars="gene.id", variable.name="SAMPLE", value.name="FPKM")}
}



phenotype.output.at.gsm.level.dt <- fread(snakemake@input[["phenotype_output_at_gsm_level_dt_filename"]])

if (FALSE) {
    phenotype.output.at.gsm.level.dt <- fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt")
}


combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.dt <-
    merge(x=combined.gexpr.FPKM.pc.only.melt.dt, y=phenotype.output.at.gsm.level.dt, by.x="SAMPLE", by.y="gsm", all.x=TRUE, all.y=FALSE) %>%
    {merge(x=., y=human.pc.gene.info.dt, by.x="gene.id", by.y="gene.id", all.x=TRUE, all.y=FALSE)} %>%
    {.[is.normal==TRUE]}

fwrite(combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.dt, snakemake@output[["combined_gexpr_FPKM_pc_only_melt_with_phenotype_normal_sample_only_dt_txt_gz_filename"]])


combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt <- combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.dt %>%
    {.[, list(FPKM.median=median(FPKM)), list(gene.id, gene.name, stage)]} %>%
    {.[, `:=`(FPKM.median.log2=log2(FPKM.median))]} %>%
    {dcast(., gene.id + gene.name ~ stage, value.var="FPKM.median.log2", fill=-Inf)} %>%
    {.[, is.maternal:=FALSE][oocyte.GV > log2(2) | oocyte.MII > log2(2), is.maternal:=TRUE]} %>%
    {.[, is.8.cell:=FALSE][`8-cell` > log2(2), is.8.cell:=TRUE]} %>%
    {.[is.maternal==TRUE, cluster:="maternal.others"][is.maternal==TRUE & pmin(oocyte.GV, oocyte.MII) > `8-cell` + 1, cluster:="maternal.decay"]} %>%
    {.[is.maternal==FALSE & is.8.cell==TRUE & `4-cell` + 1 < `8-cell`, cluster:="zygotic.genomic.activation" ]}

fwrite(combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt, snakemake@output[["combined_gexpr_FPKM_pc_only_melt_with_phenotype_normal_sample_only_median_annotated_dt_txt_gz_filename"]])
