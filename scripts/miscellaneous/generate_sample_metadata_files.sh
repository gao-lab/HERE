## 0. prepare directories

mkdir -p ./external/{DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY,DATASET_RNA_EDITING_NAME_DIRECTORY,DATASET_COLLECTION_NAME_DIRECTORY,DATASET_NAME_DIRECTORY}

## 1. Deploy GSE101571 (23 RNA-Seq samples)

raw_dataset_name="PRJNA394846_SRP112718_GSE101571.RNA.only"
dataset_name="200902-GSE101571-full"

cat ./external/NCBI.SRA.MetaData/GSE101571.txt | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{if ($2=="RNA-Seq") print $14,$1}' | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    srr=`echo $line | cut -f 2 -d ','`
    read_length_suffix=0
    if echo $srr | grep -q -E  "(SRR583734[5678]+|SRR583735[3478]+)" ; then
        read_length_suffix="125-125"
    elif echo $srr | grep -q -E  "(SRR5837349|SRR583735[012569]+|SRR583736[012]+|SRR583738[6789]+|SRR583739[0123456789]+|SRR583740[0123]+)"; then
        read_length_suffix="100-100"
    elif echo $srr | grep -q -E  "(SRR694054[45]+)" ; then
        read_length_suffix="150-150"
    else
        echo "Unsupported srr $srr" && exit 1
    fi

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done

echo ${dataset_name}-100-100 ${dataset_name}-125-125 ${dataset_name}-150-150 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-100-100
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-125-125
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-150-150

cat ./external/NCBI.SRA.MetaData/GSE101571.txt | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{if ($2=="RNA-Seq") print $14,$3}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 200 ]]; then
        type='paired-100-100'
        DATASET_NAME=${dataset_name}'-100-100'
        indexer_parameters=95
        read_length_suffix='100-100'
    elif [[ $avgspotlen == 250 ]]; then
        type='paired-125-125'
        DATASET_NAME=${dataset_name}'-125-125'
        indexer_parameters=120
        read_length_suffix='125-125'
    elif [[ $avgspotlen == 300 ]]; then
        type='paired-150-150'
        DATASET_NAME=${dataset_name}'-150-150'
        indexer_parameters=145
        read_length_suffix='150-150'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
done

## 2. Deploy GSE71318 (48 RNA-Seq samples)

raw_dataset_name="PRJNA291062_SRP061636_GSE71318"
dataset_name="200919-GSE71318-full48"

cat ./external/NCBI.SRA.MetaData/GSE71318.txt | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{if ($2=="RNA-Seq") print $14,$1,$3}' | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    srr=`echo $line | cut -f 2 -d ','`
    avgspotlen=`echo $line | cut -f 3 -d ','`
    read_length_suffix=none
    if [[ $avgspotlen == 250 ]]; then
        read_length_suffix="125-125"
    elif [[ $avgspotlen == 202 ]]; then
        read_length_suffix="101-101"
    fi
    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done

echo ${dataset_name}-101-101 ${dataset_name}-125-125 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-101-101
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-125-125


cat ./external/NCBI.SRA.MetaData/GSE71318.txt | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{if ($2=="RNA-Seq") print $14,$3}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`


    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 202 ]]; then
        type='paired-101-101'
        DATASET_NAME=${dataset_name}'-101-101'
        indexer_parameters=96
        read_length_suffix='101-101'
    elif [[ $avgspotlen == 250 ]]; then
        type='paired-125-125'
        DATASET_NAME=${dataset_name}'-125-125'
        indexer_parameters=120
        read_length_suffix='125-125'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}

done

## 3. Deploy GSE133854 (296 samples)

raw_dataset_name="PRJNA552818_SRP213116_GSE133854"
dataset_name="200924-GSE133854-all296"

cat ./external/NCBI.SRA.MetaData/GSE133854.txt | tail -n +2 | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{if ($2=="RNA-Seq") print $16,$1,$3}' | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    srr=`echo $line | cut -f 2 -d ','`
    avgspotlen=`echo $line | cut -f 3 -d ','`
    read_length_suffix=none
    if [[ $avgspotlen == 300 ]]; then
        read_length_suffix="150-150"
    elif [[ $avgspotlen == 180 ]]; then
        read_length_suffix="90-90"
    fi
    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done


echo ${dataset_name}-90-90 ${dataset_name}-150-150 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-90-90
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-150-150


cat ./external/NCBI.SRA.MetaData/GSE133854.txt | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{if ($2=="RNA-Seq") print $16,$3}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`


    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 180 ]]; then
        type='paired-90-90'
        DATASET_NAME=${dataset_name}'-90-90'
        indexer_parameters=85
        read_length_suffix='90-90'
    elif [[ $avgspotlen == 300 ]]; then
        type='paired-150-150'
        DATASET_NAME=${dataset_name}'-150-150'
        indexer_parameters=145
        read_length_suffix='150-150'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}

done


## 4. Deploy GSE136447 (508 samples)

raw_dataset_name="PRJNA562548_SRP219759_GSE136447"
dataset_name="201109-GSE136447-long508"

echo ${dataset_name}-149-149 ${dataset_name}-150-150 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-149-149
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-150-150
echo ${dataset_name}-149-149 ${dataset_name}-150-150 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-149-149
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-150-150


cat ./external/NCBI.SRA.MetaData/GSE136447.txt | tail -n +2 | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{if ($3=="RNA-Seq") print $16,$4,$1}' | grep -v ",100," | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`


    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 297 || $avgspotlen == 298 ]]; then
        type='paired-149-149'
        DATASET_NAME=${dataset_name}'-149-149'
        indexer_parameters=144
        read_length_suffix='149-149'
    elif [[ $avgspotlen == 299 ]]; then
        type='paired-150-150'
        DATASET_NAME=${dataset_name}'-150-150'
        indexer_parameters=145
        read_length_suffix='150-150'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}

    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}
    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz

done


## 5. Deploy GSE125616 (640 samples)

raw_dataset_name="PRJNA516921_SRP181926_GSE125616"
dataset_name="200911-GSE125616-all"

echo  ${dataset_name}-150-150 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-150-150

cat ./external/NCBI.SRA.MetaData/GSE125616.txt | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{if ($2=="RNA-Seq") print $15,$1}' | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    srr=`echo $line | cut -f 2 -d ','`
    read_length_suffix=150-150
    type='paired-150-150'
    DATASET_NAME=${dataset_name}'-150-150'
    indexer_parameters=145
    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
done

cat ./external/NCBI.SRA.MetaData/GSE125616.txt | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{if ($2=="RNA-Seq") print $15,$1}' | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    srr=`echo $line | cut -f 2 -d ','`
    read_length_suffix=150-150
    if  [[ ! -s /lustre4/gaog_pkuhpc/Management.dataset.repository/result/SRA-open/P01__fasterq_dump_a_single_SRA/${raw_dataset_name}/$srr/${srr}_2.fastq.gz ]]; then
        echo "$gsm / $srr not downloaded yet"
        continue
    fi
    echo $line >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}
    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done



## 6. Deploy GSE44183 (21 samples)


raw_dataset_name="PRJNA1894_SRP018525_GSE44183"
dataset_name="201217-GSE44183-earlyhumanlong21"

echo ${dataset_name}-90-90 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-90-90
echo ${dataset_name}-90-90 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-90-90

cat ./external/NCBI.SRA.MetaData/GSE44183.txt | tail -n +2 | grep 'Homo sapiens'|grep -v blood | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $15,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 180 ]]; then
        type='paired-90-90'
        DATASET_NAME=${dataset_name}'-90-90'
        indexer_parameters=85
        read_length_suffix='90-90'
    elif [[ $avgspotlen == 98 ]]; then
        continue
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done


## 7. Deploy GSE72379 (16 samples)

raw_dataset_name="PRJNA293908_SRP062850_GSE72379"
dataset_name="201101-GSE72379-full16"

echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-101-101
echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-101-101


cat ./external/NCBI.SRA.MetaData/GSE72379.txt | tail -n +2 | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $16,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 202 ]]; then
        type='paired-101-101'
        DATASET_NAME=${dataset_name}'-101-101'
        indexer_parameters=96
        read_length_suffix='101-101'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done


## 8. Deploy GSE36552 (124 samples)

raw_dataset_name="PRJNA153427_SRP011546_GSE36552"
dataset_name="201104-GSE36552-full124"

echo ${dataset_name}-100 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-100
echo ${dataset_name}-100 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-100


cat ./external/NCBI.SRA.MetaData/GSE36552.txt | tail -n +2 | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $16,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 100 ]]; then
        type='single-100'
        DATASET_NAME=${dataset_name}'-100'
        indexer_parameters=95
        read_length_suffix='100'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r.fastq.gz
done

## We found that there is one sample () that contains invalid reads
## We decided to manually remove invalid reads from this sample.

cd external/RNA-Seq-with-Run/201104-GSE36552-full124-100/GSM922196/SRR491011/RNA
mv r.fastq.gz r.raw.fastq.gz
## only one record: SRR491011.28354271
## See https://www.biostars.org/p/66996/ for more details of this filtering
zcat r.raw.fastq.gz | awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) == length(qseq)) {print header, seq, qheader, qseq}}' | gzip > r.filtered.fastq.gz
ln -s r.filtered.fastq.gz r.fastq.gz
cd -


## 9. Deploy GSE95477 (20 samples)

raw_dataset_name="PRJNA377237_SRP100829_GSE95477"
dataset_name="201101-GSE95477-full20"

echo ${dataset_name}-100-100 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-100-100

echo ${dataset_name}-100-100 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-100-100



cat ./external/NCBI.SRA.MetaData/GSE95477.txt | tail -n +2 | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $16,$4,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 200 ]]; then
        type='paired-100-100'
        DATASET_NAME=${dataset_name}'-100-100'
        indexer_parameters=95
        read_length_suffix='100-100'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done



## 10. Deploy GSE65481 (22 samples)

raw_dataset_name="PRJNA274140_SRP053004_GSE65481"
dataset_name="201031-GSE65481-full22"

echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-101-101
echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-101-101


cat ./external/NCBI.SRA.MetaData/GSE65481.txt | tail -n +2 | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $15,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 202 ]]; then
        type='paired-101-101'
        DATASET_NAME=${dataset_name}'-101-101'
        indexer_parameters=96
        read_length_suffix='101-101'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done

## 11. Deploy GSE130289 (139 samples)

raw_dataset_name="PRJNA534673_SRP193790_GSE130289"
dataset_name="201031-GSE130289-full139"

echo ${dataset_name}-101-101 ${dataset_name}-150-150 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-101-101
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-150-150
echo ${dataset_name}-101-101 ${dataset_name}-150-150 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-101-101
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-150-150


cat ./external/NCBI.SRA.MetaData/GSE130289.txt | tail -n +2 | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $15,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 202 ]]; then
        type='paired-101-101'
        DATASET_NAME=${dataset_name}'-101-101'
        indexer_parameters=96
        read_length_suffix='101-101'
    elif [[ $avgspotlen == 300 ]]; then
        type='paired-150-150'
        DATASET_NAME=${dataset_name}'-150-150'
        indexer_parameters=145
        read_length_suffix='150-150'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done


## 12. Deploy GSE100118 (92 samples)


raw_dataset_name="PRJNA392971_SRP113531_GSE100118"
dataset_name="201101-GSE100118-full92"

echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-101-101
echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-101-101


cat ./external/NCBI.SRA.MetaData/GSE100118.txt | tail -n +2 | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $16,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 202 ]]; then
        type='paired-101-101'
        DATASET_NAME=${dataset_name}'-101-101'
        indexer_parameters=96
        read_length_suffix='101-101'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done



## 13. Deploy GSE49828 (3 samples)


raw_dataset_name="PRJNA215030_SRP028804_GSE49828.RNA.Seq.only"
dataset_name="201104-GSE49828-RNASeqonly3"

echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-101-101
echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-101-101


cat ./external/NCBI.SRA.MetaData/GSE49828.txt | tail -n +2 | grep "RNA-Seq" | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $15,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 202 ]]; then
        type='paired-101-101'
        DATASET_NAME=${dataset_name}'-101-101'
        indexer_parameters=96
        read_length_suffix='101-101'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done



## 14. Deploy GSE64417 (21 hESC samples)


raw_dataset_name="PRJNA270993_SRP051472_GSE64417"
dataset_name="201218-GSE64417-hESConly21"

echo ${dataset_name}-100-100 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-100-100
echo ${dataset_name}-100-100 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-100-100



cat ./external/NCBI.SRA.MetaData/GSE64417.txt | tail -n +2 | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $14,$3,$1}' | grep -P  "GSM15710(5|6|7[012345]+)" | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 200 ]]; then
        type='paired-100-100'
        DATASET_NAME=${dataset_name}'-100-100'
        indexer_parameters=95
        read_length_suffix='100-100'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done


## 15. Deploy GSE62772 (18 hESC samples)

raw_dataset_name="PRJNA265099_SRP049340_GSE62772.hESC.only"
dataset_name="201102-GSE62772-hESC18"

echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-101-101
echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-101-101


cat ./external/NCBI.SRA.MetaData/GSE62772.txt | tail -n +2 | grep "reference hESC line" | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $14,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 202 ]]; then
        type='paired-101-101'
        DATASET_NAME=${dataset_name}'-101-101'
        indexer_parameters=96
        read_length_suffix='101-101'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done



## 16. Deploy GSE126488 (40 samples)

raw_dataset_name="PRJNA522065_SRP185781_GSE126488"
dataset_name="201103-GSE126488-full40"

echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-101-101
echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-101-101


cat ./external/NCBI.SRA.MetaData/GSE126488.txt | tail -n +2 | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $16,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 202 ]]; then
        type='paired-101-101'
        DATASET_NAME=${dataset_name}'-101-101'
        indexer_parameters=96
        read_length_suffix='101-101'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done



## 17. Deploy GSE73211 (30 hESC samples)

raw_dataset_name="PRJNA296379_SRP063867_GSE73211"
dataset_name="201102-GSE73211-ESC35"

echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-101-101
echo ${dataset_name}-101-101 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-101-101


cat ./external/NCBI.SRA.MetaData/GSE73211.txt | tail -n +2 | grep "ESC" | grep 'RNA-Seq' | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $15,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 202 ]]; then
        type='paired-101-101'
        DATASET_NAME=${dataset_name}'-101-101'
        indexer_parameters=96
        read_length_suffix='101-101'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done


## 18. Deploy GSE119324 (10 samples)

raw_dataset_name="PRJNA488795_SRP159258_GSE119324"
dataset_name="201104-GSE119324-full10"

echo ${dataset_name}-100-100 | tr ' ' '\n' > ./external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo TYPE,DATASET_NAME,SAMPLE_NAME,INDEXER_PARAMETERS  > ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${dataset_name}-100-100
echo ${dataset_name}-100-100 | tr ' ' '\n' > ./external/DATASET_COLLECTION_NAME_DIRECTORY/${dataset_name}
echo SAMPLE_NAME,RUN_NAME  > ./external/DATASET_NAME_DIRECTORY/${dataset_name}-100-100


cat ./external/NCBI.SRA.MetaData/GSE119324.txt | tail -n +2 | awk -v FPAT='[^,]*|"[^"]+"' -v OFS="," '{ print $15,$3,$1}' | sort | uniq | while read line
do
    gsm=`echo $line | cut -f 1 -d ','`
    avgspotlen=`echo $line | cut -f 2 -d ','`
    srr=`echo $line | cut -f 3 -d ','`

    type=none
    DATASET_NAME=none
    indexer_parameters=none
    read_length_suffix=none
    if [[ $avgspotlen == 200 ]]; then
        type='paired-100-100'
        DATASET_NAME=${dataset_name}'-100-100'
        indexer_parameters=95
        read_length_suffix='100-100'
    fi

    echo $type,$DATASET_NAME,$gsm,$indexer_parameters >> ./external/DATASET_RNA_EDITING_NAME_DIRECTORY/${DATASET_NAME}
    echo $gsm,$srr >> ./external/DATASET_NAME_DIRECTORY/${dataset_name}-${read_length_suffix}

    rm -fr ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/*
    mkdir -p ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/
##  ln -s YOUR-PATH-with-${srr}-TO-r1.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r1.fastq.gz
##  ln -s YOUR-PATH-with-${srr}-TO-r2.fastq.gz ./external/RNA-Seq-with-Run/${dataset_name}-${read_length_suffix}/$gsm/$srr/RNA/r2.fastq.gz
done
