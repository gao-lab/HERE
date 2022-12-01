import pandas
import re
import os



def __get_parameters_collection_for_RNA_editing_calling(dataset_collection_name, prefix_for_dataset_name_directory="external/DATASET_RNA_EDITING_NAME_DIRECTORY/", prefix_for_dataset_collection_name_directory="external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/"):
    dataset_name_collection = __get_dataset_name_collection(dataset_collection_name, prefix_for_dataset_collection_name_directory=prefix_for_dataset_collection_name_directory)
    final_collection = []
    for temp_dataset_name in dataset_name_collection:
        temp_dataset_DataFrame = pandas.read_csv( prefix_for_dataset_name_directory + temp_dataset_name, sep=",")
        temp_row_collection = [ temp_row for temp_rowindex, temp_row in temp_dataset_DataFrame.iterrows()]
        final_collection = final_collection + temp_row_collection
    return final_collection



def __get_dataset_name_collection(dataset_collection_name, prefix_for_dataset_collection_name_directory="external/DATASET_COLLECTION_NAME_DIRECTORY/"):
    dataset_name_collection = []
    with open(prefix_for_dataset_collection_name_directory + dataset_collection_name, "r") as temp_fileobject:
        dataset_name_collection = [temp_filename.strip() for temp_filename in temp_fileobject.readlines() if temp_filename != '\n']
    return dataset_name_collection


def __get_dataset_and_all_its_samples_name_collection(dataset_collection_name, prefix_for_dataset_name_directory="external/DATASET_NAME_DIRECTORY/", prefix_for_dataset_collection_name_directory="external/DATASET_COLLECTION_NAME_DIRECTORY/"):
    dataset_name_collection = __get_dataset_name_collection(dataset_collection_name, prefix_for_dataset_collection_name_directory=prefix_for_dataset_collection_name_directory)
    final_collection = []
    for temp_dataset_name in dataset_name_collection:
        temp_dataset_DataFrame = pandas.read_csv( prefix_for_dataset_name_directory + temp_dataset_name, sep=",")
        temp_sample_name_collection = temp_dataset_DataFrame.loc[:, "SAMPLE_NAME"]
        temp_dataset_and_all_its_samples_name_collection = [ temp_dataset_name + "/" + temp_sample_name for temp_sample_name in temp_sample_name_collection]
        final_collection = final_collection + temp_dataset_and_all_its_samples_name_collection
    return final_collection


def __get_dataset_and_all_its_samples_and_all_their_runs_name_collection(dataset_collection_name, prefix_for_dataset_name_directory="external/DATASET_NAME_DIRECTORY/", prefix_for_dataset_collection_name_directory="external/DATASET_COLLECTION_NAME_DIRECTORY/"):
    dataset_name_collection = __get_dataset_name_collection(dataset_collection_name, prefix_for_dataset_collection_name_directory=prefix_for_dataset_collection_name_directory)
    final_collection = []
    for temp_dataset_name in dataset_name_collection:
        temp_dataset_DataFrame = pandas.read_csv( prefix_for_dataset_name_directory + temp_dataset_name, sep=",")
        temp_sample_name_collection = temp_dataset_DataFrame.loc[:, "SAMPLE_NAME"]
        temp_run_name_collection = temp_dataset_DataFrame.loc[:, "RUN_NAME"]
        temp_dataset_and_all_its_samples_and_all_its_runs_name_collection = [ temp_dataset_name + "/" + temp_sample_name + "/" + temp_run_name for temp_sample_name, temp_run_name in zip(temp_sample_name_collection, temp_run_name_collection)]
        final_collection = final_collection + temp_dataset_and_all_its_samples_and_all_its_runs_name_collection
    return final_collection


rule S90_1__check_recal_bam_stat:
    input:
        flag_S15_1__step02__part05="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part05__locally_realign_reads"
    output:
        flag_S90_1=touch("result/S90_1__check_recal_bam_stat/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished"),
        stats="result/S90_1__check_recal_bam_stat/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/stats",
        coverage="result/S90_1__check_recal_bam_stat/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/coverage"
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/"
    threads:
        1
    shell:
        """
        samtools stats "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam" > "{output.stats}"
        samtools coverage "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam" > "{output.coverage}"
        """

rule B90_1__check_recal_bam_stat:
    input:
        lambda wildcards: ["result/S90_1__check_recal_bam_stat/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/" + str(temp_row.INDEXER_PARAMETERS) + "/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_B90_1=touch("result/B90_1__check_recal_bam_stat/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished")



rule S92_1__summarize_sample_stat:
    input:
        flag_S90_1="result/S90_1__check_recal_bam_stat/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished",
        stats="result/S90_1__check_recal_bam_stat/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/stats",
        coverage="result/S90_1__check_recal_bam_stat/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/coverage",
        flag_S91_1="result/S91_1__merge_trim_summary/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/finished",
        merged_trim_summary="result/S91_1__merge_trim_summary/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/merged.trim.summary"
    output:
        flag_S92_1=touch("result/S92_1__summarize_sample_stat/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished"),
        summary_dt="result/S92_1__summarize_sample_stat/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/summary.dt.gz"
    threads:
        1
    script:
        "./scripts.for.report.ver2/basic.summary/run_summarize_bam_stats.R"


rule B92_1__summarize_sample_stat:
    input:
        flags_collection = lambda wildcards: ["result/S92_1__summarize_sample_stat/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/" + str(temp_row.INDEXER_PARAMETERS) + "/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)],
        summary_dts_collection = lambda wildcards: ["result/S92_1__summarize_sample_stat/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/" + str(temp_row.INDEXER_PARAMETERS) + "/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/summary.dt.gz" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_B92_1=touch("result/B92_1__summarize_sample_stat/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished"),
        combined_summary_dt="result/B92_1__summarize_sample_stat/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/combined.summary.dt.gz"
    script:
        "./scripts.for.report.ver2/basic.summary/REDUCE.run_summarize_bam_stats.R"

