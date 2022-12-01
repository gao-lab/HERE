import pandas
import re
import os


rule S80_1__compute_AEI:
    input:
        flag_S15_1__step02__part05="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part05__locally_realign_reads"
    output:
        flag_S80_1__compute_AEI=touch("result/S80_1__compute_AEI/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished"),
        Editing_index_csv="result/S80_1__compute_AEI/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/EditingIndex.csv"
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        AEI_result_directory="result/S80_1__compute_AEI/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/"
    conda:
        "RNAEditingIndexer"
    threads:
        config['threads_computing_AEI']
    shell:
        """
        rm -fr "{params.AEI_result_directory}"/flags
        rm -f "{params.AEI_result_directory}"/temp.input/linked.bam
        mkdir -p "{params.AEI_result_directory}"/temp.input
        ln -s -r "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.realn.bam" "{params.AEI_result_directory}"/temp.input/linked.bam
        rm -fr "{params.AEI_result_directory}"/temp.output/*
        mkdir -p "{params.AEI_result_directory}"/temp.output
        rm -fr "{params.AEI_result_directory}"/temp.output.summary/*
        mkdir -p "{params.AEI_result_directory}"/temp.output.summary
        rm -f "{params.AEI_result_directory}"/RNAEditingIndex
        ln -s -r ./tools/RNAEditingIndexer/RNAEditingIndex "{params.AEI_result_directory}"/RNAEditingIndex

        pairedflag="--paired_end"
        type_end=`echo {wildcards.TYPE} | cut -f 1 -d '-'`
        if [[ $type_end == 'single' ]]; then
            pairedflag=""
        fi

        cd "{params.AEI_result_directory}"
        current_absolute_path=`pwd`
        ./RNAEditingIndex $pairedflag -d $current_absolute_path/temp.input -f linked.bam --output_dir $current_absolute_path/temp.output  --output_dir_summery $current_absolute_path/temp.output.summary --genome hg38 --ts 1 --tsd {threads}
        ln -s -r temp.output.summary/EditingIndex.csv EditingIndex.csv
        cd -
        """

def __get_dataset_name_collection(dataset_collection_name, prefix_for_dataset_collection_name_directory="external/DATASET_COLLECTION_NAME_DIRECTORY/"):
    dataset_name_collection = []
    with open(prefix_for_dataset_collection_name_directory + dataset_collection_name, "r") as temp_fileobject:
        dataset_name_collection = [temp_filename.strip() for temp_filename in temp_fileobject.readlines() if temp_filename != '\n']
    return dataset_name_collection
        
        
def __get_parameters_collection_for_RNA_editing_calling(dataset_collection_name, prefix_for_dataset_name_directory="external/DATASET_RNA_EDITING_NAME_DIRECTORY/", prefix_for_dataset_collection_name_directory="external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/"):
    dataset_name_collection = __get_dataset_name_collection(dataset_collection_name, prefix_for_dataset_collection_name_directory=prefix_for_dataset_collection_name_directory)
    final_collection = []
    for temp_dataset_name in dataset_name_collection:
        temp_dataset_DataFrame = pandas.read_csv( prefix_for_dataset_name_directory + temp_dataset_name, sep=",")
        temp_row_collection = [ temp_row for temp_rowindex, temp_row in temp_dataset_DataFrame.iterrows()]
        final_collection = final_collection + temp_row_collection
    return final_collection


rule BS80_1__compute_AEI:
    input:
        lambda wildcards: ["result/S80_1__compute_AEI/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/" + str(temp_row.INDEXER_PARAMETERS) + "/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)],
        Editing_file_collection = lambda wildcards: list(set(["result/S80_1__compute_AEI/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/" + str(temp_row.INDEXER_PARAMETERS) + "/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/EditingIndex.csv" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]))
    output:
        flag_B15_1__step10=touch("result/BS80_1__compute_AEI/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished"),
        combined_Editing_file="result/BS80_1__compute_AEI/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/combined.Editing.file.dt.gz"
    script:
        "./scripts/BS80_1__compute_AEI/run.R"
