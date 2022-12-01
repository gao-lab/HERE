library("data.table")
library("plyr")
library("readxl")


input.phenotype_collection_filename <- snakemake@input[["phenotype_collection_filename"]]
input.manuscript_table_for_processing_xlsx_filename  <- snakemake@input[["manuscript_table_for_processing_xlsx_filename"]]

output.phenotype_output_dt_filename <- snakemake@output[["phenotype_output_dt_filename"]]
output.phenotype_output_at_gsm_level_dt_filename <- snakemake@output[["phenotype_output_at_gsm_level_dt_filename"]]




phenotype.collection.dt <- fread(input.phenotype_collection_filename)
manuscript.table.for.processing.dt <- data.table(read_excel(input.manuscript_table_for_processing_xlsx_filename, sheet="phenotype"))

phenotype.output.dt <- rbindlist(alply(phenotype.collection.dt, .margins=1, .fun=function(temp.row.dt){
    temp.gse <- temp.row.dt[1, STUDY]
    temp.processing.row.dt <- manuscript.table.for.processing.dt[GSE == temp.gse]
    temp.gse.and.sra.dt <- fread(paste(sep="", "./external/NCBI.SRA.MetaData/", temp.row.dt[1, PHENOTYPE_FILENAME]))[, gsm:=`Sample Name`]


    temp.gse.and.gsm.dt <- fread(paste(sep="", "result/S20_1__extract_GSE_table_by_GEOmetadb/", temp.gse, "/gse.and.gsm.dt.txt"))
    temp.gse.and.sra.and.gsm.dt <- copy(temp.gse.and.sra.dt)
    if (nrow(temp.gse.and.gsm.dt) > 0){
        temp.gse.and.sra.and.gsm.dt <- merge(x=temp.gse.and.sra.dt, y=temp.gse.and.gsm.dt, by.x="gsm", by.y="gsm", all.x=TRUE, all.y=FALSE)
    }

    temp.output.dt <- temp.gse.and.sra.and.gsm.dt[eval(parse(text=temp.processing.row.dt[1, filter]))][, list(gse=temp.gse, gsm, srr=Run, srr.avgspotlen=AvgSpotLen, srr.bytes=Bytes, srr.bases=Bases, srr.avgreadcount=Bases/AvgSpotLen, stage=eval(parse(text=temp.processing.row.dt[1, stage])), is.normal=eval(parse(text=temp.processing.row.dt[1, is.normal])), disease=eval(parse(text=temp.processing.row.dt[1, disease])), treatment=eval(parse(text=temp.processing.row.dt[1, treatment])), maternal.age=eval(parse(text=temp.processing.row.dt[1, maternal.age])), developmental.day=eval(parse(text=temp.processing.row.dt[1, developmental.day])), cell.line=eval(parse(text=temp.processing.row.dt[1, cell.line])) )]


    return(temp.output.dt)
}, .progress='text'))

phenotype.output.at.gsm.level.dt <- phenotype.output.dt[, list(srr.count=.N, srr.mean.avgspotlen=mean(srr.avgspotlen), srr.total.bytes=sum(srr.bytes), srr.total.bases=sum(srr.bases), srr.total.avgreadcount=sum(srr.avgreadcount)), list(gse, gsm, stage, is.normal, disease, treatment, maternal.age, developmental.day, cell.line)]

fwrite(phenotype.output.dt, output.phenotype_output_dt_filename)
fwrite(phenotype.output.at.gsm.level.dt, output.phenotype_output_at_gsm_level_dt_filename)
