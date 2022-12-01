#!/bin/bash

tar -czvf ./zenodo-archives-CB.rev1/UCSC.Tables.tar.gz \
    external/UCSC.Table.Browser.knownGene.GENCODE/32/knownGene \
    external/UCSC.Table.Browser.dbSNP/151/flagged.cDNA.only/dbSNP.bed \
    external/UCSC.Table.Browser.repeatmasker/repFamily.Alu/repeatmasker.bed \
    external/UCSC.Table.Browser.repeatmasker/repFamily.Simple_repeat/repeatmasker.bed \
    external/UCSC.Table.Browser.repeatmasker/repFamily.not.Alu/repeatmasker.bed
