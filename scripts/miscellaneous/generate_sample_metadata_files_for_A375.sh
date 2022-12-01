## 0. prepare directories

mkdir -p ./external/{DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY,DATASET_RNA_EDITING_NAME_DIRECTORY,DATASET_COLLECTION_NAME_DIRECTORY,DATASET_NAME_DIRECTORY}

## 1. Deploy GSE144296 A375 (DNA)

#+BEGIN_SRC sh
raw_dataset_name="GSE144296.A375"
dataset_name="210203-GSE144296.A375-DNA"


echo ${dataset_name}-37-37 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-37-37

echo ${dataset_name}-37-37 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-37-37

cat ./external/NCBI.SRA.MetaData/GSE144296.txt | tail -n +2 | grep -P ",A375," | grep -P ",OTHER," | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $18,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters="DNA"
    read_length_suffix=none
    if [[ $avgspotlen == 74 ]]; then
        type='paired-37-37'
        DATASET_NAME=${dataset_name}'-37-37'
        read_length_suffix='37-37'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done

#+END_SRC


## 2. Deploy GSE144296 A375 (RNA)

#+BEGIN_SRC sh
raw_dataset_name="GSE144296.A375"
dataset_name="210203-GSE144296.A375-RNA"


echo ${dataset_name}-37-37 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-37-37

echo ${dataset_name}-37-37 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-37-37

cat ./external/NCBI.SRA.MetaData/GSE144296.txt | tail -n +2 | grep -P ",A375," | grep -P ",RNA-Seq," | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $18,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 74 ]]; then
        type='paired-37-37'
        DATASET_NAME=${dataset_name}'-37-37'
        indexer_parameters=32
        read_length_suffix='37-37'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done

#+END_SRC

