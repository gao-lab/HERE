library("GEOmetadb")
library("data.table")

input.GEOmetadb_sqlite_filename <- snakemake@input[["GEOmetadb_sqlite_filename"]]
wildcards.gse <- snakemake@wildcards[["GSE"]]
output.gse_and_gsm_dt_filename <- snakemake@output[["gse_and_gsm_dt_filename"]]

GEO.connection <- dbConnect(SQLite(), input.GEOmetadb_sqlite_filename)

gse.and.gsm.dt <- data.table(dbGetQuery(GEO.connection, paste(sep="", "SELECT ", ' gse.gse as gse,
gse.title as "gse.title",
gse.status as "gse.status",
gse.submission_date as "gse.submission.date",
gse.last_update_date as "gse.last.update.date",
gse.pubmed_id as "gse.pubmed.id",
gse.summary as "gse.summary",
gse.type as "gse.type",
gse.contributor as "gse.contributor",
gse.overall_design as "gse.overall.design",
gse.repeats as "gse.repeats",
gse.repeats_sample_list as "gse.repeats.sample.list",
gse.variable as "gse.variable",
gse.variable_description as "gse.variable.descrciption",
gse.contact as "gse.contact",
gse.supplementary_file as "gse.supplementary.file",
gse_gsm.gsm as gsm,
gsm.title as "gsm.title",
gsm.gpl as "gsm.gpl",
gsm.status as "gsm.status",
gsm.submission_date as "gsm.submission.date",
gsm.last_update_date as "gsm.last.update.date",
gsm.type as "gsm.type",
gsm.source_name_ch1 as "gsm.source.name.ch1",
gsm.organism_ch1 as "gsm.organism.ch1",
gsm.characteristics_ch1 as "gsm.characteristics.ch1",
gsm.molecule_ch1 as "gsm.molecule.ch1",
gsm.label_ch1 as "gsm.label.ch1",
gsm.treatment_protocol_ch1 as "gsm.treatment.protocol.ch1",
gsm.extract_protocol_ch1 as "gsm.extract.protocol.ch1",
gsm.label_protocol_ch1 as "gsm.label.protocol.ch1",
gsm.source_name_ch2 as "gsm.source.name.ch2",
gsm.organism_ch2 as "gsm.organism.ch2",
gsm.characteristics_ch2 as "gsm.characteristics.ch2",
gsm.molecule_ch2 as "gsm.molecule.ch2",
gsm.label_ch2 as "gsm.label.ch2",
gsm.treatment_protocol_ch2 as "gsm.treatment.protocol.ch2",
gsm.extract_protocol_ch2 as "gsm.extract.protocol.ch2",
gsm.label_protocol_ch2 as "gsm.label.protocol.ch2",
gsm.hyb_protocol as "gsm.hyb.protocol",
gsm.description as "gsm.description",
gsm.data_processing as "gsm.data.processing",
gsm.contact as "gsm.contact",
gsm.supplementary_file as "gsm.supplementary.file",
gsm.data_row_count as "gsm.data.row.count",
gsm.channel_count as "gsm.channel.count"
  FROM gse INNER JOIN gse_gsm ON gse.gse = gse_gsm.gse INNER JOIN gsm ON gse_gsm.gsm = gsm.gsm WHERE gse.gse == "', wildcards.gse, '"')))


fwrite(gse.and.gsm.dt, output.gse_and_gsm_dt_filename)

dbDisconnect(GEO.connection)
