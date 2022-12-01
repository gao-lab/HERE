library(data.table)
library(magrittr)

args <- commandArgs(trailingOnly=TRUE)

input.filename <- args[1]
taxonomy.filename <- args[2]
output.filename <- args[3]

maf.table.dt <- fread(input.filename, header=FALSE, col.names=c('build', 'transcript.id', 'alignment'))

taxonomy.build.to.id.mapping.vector <- fread(taxonomy.filename, header=FALSE, col.names=c('build', 'id', 'group')) %>% {{temp.vector <- .[, id]; names(temp.vector) <- .[, build]; temp.vector}}

maf.table.dt[, id:=taxonomy.build.to.id.mapping.vector[build]]

fwrite(maf.table.dt[, list(transcript.id, id, alignment)], output.filename, col.names=FALSE, sep='\t')
