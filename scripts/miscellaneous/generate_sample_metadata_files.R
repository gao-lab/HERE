library("data.table")
library("glue")
library("foreach")
library("iterators")

dir.create("external/DATASET_NEW_NAME_DIRECTORY", recursive=TRUE)
dir.create("external/DATASET_NEW_COLLECTION_NAME_DIRECTORY", recursive=TRUE)

{    
    "./scripts/miscellaneous/info.csv" -> .;
    fread(.) -> .
    .[, ID:=GSE] -> .;
    ## load NCBI SRA MetaData
    rbindlist(foreach(TEMP.ROW.DT=iter(., by="row")) %do% {
        fread(glue("external/NCBI.SRA.MetaData/{TEMP.ROW.DT[1, GSE]}.txt")) -> TEMP.DT;
        data.table(TEMP.DT, TEMP.ROW.DT) -> TEMP.RESULT.DT
        TEMP.RESULT.DT
    }, use.names=TRUE, fill=TRUE) -> .;
    ##
    ## Filter samples
    ## 1. RNA-Seq only
    .[`Assay Type`=='RNA-Seq'] -> .;
    ## 2. For datasets other than GSE36552, no short AvgSpotLen 
    .[ (GSE=="GSE36552") | ((LibraryLayout == "PAIRED") & AvgSpotLen >= 150)] -> .;
    ## 3. Pick relevant embryos for datasets with mixed samples
    .[( (GSE!="GSE44183") | (Organism == "Homo sapiens") ) == TRUE] -> .;
    .[( (GSE!="GSE64417") | (grepl("GSM15710(5|6|7[012345]+)", `Sample Name`)) ) == TRUE] -> .;
    .[( (GSE!="GSE62772") | (source_name == "reference hESC line") ) == TRUE] -> .;
    .[( (GSE!="GSE73211") | (Cell_type == "ESC") ) == TRUE] -> .;    
    .[((GSE=="GSE136447") & (AvgSpotLen == 100)) == FALSE] -> .;
    .[((GSE=="GSE44183") & (
        (Organism != "Homo sapiens") |
        (AvgSpotLen == 98))
    )==FALSE] -> .;
    ## add indexer parameters and DATASET_NAME
    .[GSE=="GSE36552", BWA_INDEXER_PARAMETERS:=95] -> .;
    .[GSE=="GSE36552", STAR_INDEXER_PARAMETERS:=99] -> .;
    .[GSE=="GSE36552", DATASET_NAME_SUFFIX:="100"] -> .;
    .[GSE=="GSE36552", TYPE:=paste(sep="", "single-", DATASET_NAME_SUFFIX)] -> .;
    .[GSE!="GSE36552", BWA_INDEXER_PARAMETERS:=c("180"=85, "200"=95, "202"=96, "250"=120, "297"=144, "298"=144, "299"=145, "300"=145)[as.character(AvgSpotLen)]] -> .;
    .[GSE!="GSE36552", STAR_INDEXER_PARAMETERS:=c("180"=89, "200"=99, "202"=100, "250"=124, "297"=148, "298"=148, "299"=149, "300"=149)[as.character(AvgSpotLen)]] -> .;
    .[GSE!="GSE36552", DATASET_NAME_SUFFIX:=c("180"="90-90", "200"="100-100", "202"="101-101", "250"="125-125", "297"="149-149", "298"="149-149", "299"="150-150", "300"="150-150")[as.character(AvgSpotLen)]]
    .[GSE!="GSE36552", TYPE:=paste(sep="", "paired-", DATASET_NAME_SUFFIX)] -> .;
    .[, DATASET_NAME:=paste(sep="", ROUGH_DATASET_NAME, "-", DATASET_NAME_SUFFIX)] -> .;
    ## keep necessary columns
    .[, list(GSE, TYPE, DATASET_NAME, SAMPLE_NAME=`Sample Name`, RUN_NAME=Run, BWA_INDEXER_PARAMETERS, STAR_INDEXER_PARAMETERS)][, Group:=DATASET_NAME] -> .;
    ##
    ## write meta to file
    .[, fwrite(.SD[, list(TYPE, DATASET_NAME, SAMPLE_NAME, BWA_INDEXER_PARAMETERS, STAR_INDEXER_PARAMETERS, RUN_NAME)], file=paste(sep="", "./external/DATASET_NEW_NAME_DIRECTORY/", .SD[1, DATASET_NAME]), row.names=FALSE, col.names=TRUE, sep=",", quote=FALSE), list(Group)] -> not.used.variable;
    ## write dataset
    writeLines(.[, unique(sort(DATASET_NAME))], "./external/DATASET_NEW_COLLECTION_NAME_DIRECTORY/210215-sixth-dataset")
    ##
}
