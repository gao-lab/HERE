import pandas
import re
import os

def __get_RUNs_by_DATASET_and_SAMPLE(temp_dataset_name, temp_sample_name, prefix_for_dataset_name_directory="external/DATASET_NAME_DIRECTORY/"):
    temp_dataset_DataFrame = pandas.read_csv( prefix_for_dataset_name_directory + temp_dataset_name, sep=",")
    ## print(temp_dataset_DataFrame)
    temp_run_name_for_the_sample_collection = temp_dataset_DataFrame.loc[temp_dataset_DataFrame["SAMPLE_NAME"] == temp_sample_name, "RUN_NAME"]
    ## print(temp_run_name_for_the_sample_collection)
    if (len(temp_run_name_for_the_sample_collection) == 0):
        raise ValueError("DATASET: " + temp_dataset_name + " - SAMPLE: " + temp_sample_name + " has no runs available" )
    temp_run_names_list = [ temp_run_name for temp_run_name in temp_run_name_for_the_sample_collection]
    return temp_run_names_list

def __get_trim_summaries_collection_by_DATASET_and_SAMPLE(wildcards):
    trim_summary_basenames = None
    if wildcards.TYPE.startswith("single"):
        trim_summary_basenames=["r.fastq.gz_trimming_report.txt"]
    elif wildcards.TYPE.startswith("paired"):
        trim_summary_basenames=["r1.fastq.gz_trimming_report.txt", "r2.fastq.gz_trimming_report.txt"]
    else:
        raise ValueError("invalid " + wildcards.TYPE)
    all_runs_collection=__get_RUNs_by_DATASET_and_SAMPLE(wildcards.DATASET, wildcards.SAMPLE)
    ## print(all_runs_collection)
    trim_summaries_collection = ["result/S01_2__trim_RNA_Seq/" + wildcards.TYPE + "/" + wildcards.DATASET + "/" + wildcards.SAMPLE + "/" + temp_run + "/" + wildcards.FASTP_PARAMETERS + "/" + temp_basename for temp_run in all_runs_collection for temp_basename in trim_summary_basenames]
    return trim_summaries_collection

rule S91_1__merge_trim_summary:
    input:
        flags_collection_S15_1__step01 = lambda wildcards: ["result/S01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/" + temp_run + "/{FASTP_PARAMETERS}/finished.step01__trim_RNA_Seq" for temp_run in __get_RUNs_by_DATASET_and_SAMPLE(wildcards.DATASET, wildcards.SAMPLE)]
    output:
        flag_S91_1=touch("result/S91_1__merge_trim_summary/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/finished"),
        merged_trim_summary="result/S91_1__merge_trim_summary/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/merged.trim.summary"
    params:
        trim_summaries_collection=lambda wildcards: __get_trim_summaries_collection_by_DATASET_and_SAMPLE(wildcards)
    threads:
        1
    shell:
        """
        for temp in {params.trim_summaries_collection:q}
        do
            echo "BEGINNING_OF_SUMMARY_FILE: "$temp
            cat $temp
            echo "END_OF_SUMMARY_FILE: "$temp
        done > {output.merged_trim_summary}
        """

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


rule B91_1__merge_trim_summary:
    input:
        lambda wildcards: ["result/S91_1__merge_trim_summary/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/{FASTP_PARAMETERS}/finished" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_B91_1=touch("result/B91_1__merge_trim_summary/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/finished")
