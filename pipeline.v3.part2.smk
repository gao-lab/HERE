import pandas
import re
import os

for temp_key in ['threads_indexing', 'threads_trimming', 'threads_aligning', 'threads_calling_variants', 'threads_merging_bams', 'threads_calling_expression', 'threads_merging_vcfs']:
    if temp_key not in config.keys():
        config[temp_key] = 20

for temp_key in ['threads_auxiliary_processing', 'threads_predicting_secondary_structure', 'threads_calling_variants_minor', 'threads_annotating']:
    if temp_key not in config.keys():
        config[temp_key] = 1


for temp_key in ['threads_auxiliary_processing_parallel']:
    if temp_key not in config.keys():
        config[temp_key] = 6


for temp_key in ['python2_exec']:
    if temp_key not in config.keys():
        config[temp_key] = 'none'

if 'default_Xmx' not in config.keys():
    config['default_Xmx']='-Xmx150G'

for temp_key in ['convertCoordinates_Xmx', 'GATK_RealignerTargetCreator_Xmx', 'GATK_IndelRealigner_Xmx', 'GATK_BaseRecalibrator_Xmx', 'GATK_PrintReads_Xmx', 'GATK_UnifiedGenotyper_Xmx', 'GATK_VariantFiltration_Xmx', 'picard_Xmx', 'snpEff_Xmx']:
    if temp_key not in config.keys():
        config[temp_key] = config['default_Xmx']


def __get_parameters_collection_for_RNA_editing_calling(dataset_collection_name, prefix_for_dataset_name_directory="external/DATASET_RNA_EDITING_NAME_DIRECTORY/", prefix_for_dataset_collection_name_directory="external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/"):
    dataset_name_collection = __get_dataset_name_collection(dataset_collection_name, prefix_for_dataset_collection_name_directory=prefix_for_dataset_collection_name_directory)
    final_collection = []
    for temp_dataset_name in dataset_name_collection:
        temp_dataset_DataFrame = pandas.read_csv( prefix_for_dataset_name_directory + temp_dataset_name, sep=",")
        temp_row_collection = [ temp_row for temp_rowindex, temp_row in temp_dataset_DataFrame.iterrows()]
        final_collection = final_collection + temp_row_collection
    return final_collection

def __write_temp_vcf_path_collection_to_file_and_return_filename(dataset_collection_name, temp_vcf_path_collection_filename, wildcards, prefix_for_dataset_name_directory="external/DATASET_RNA_EDITING_NAME_DIRECTORY/", prefix_for_dataset_collection_name_directory="external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/"):
    temp_directory = os.path.dirname(temp_vcf_path_collection_filename)
    os.makedirs(temp_directory, exist_ok=True)
    temp_vcf_path_collection = set(["result/S15_1__get_sample_RNA_editing_sites_v3/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/"  + str(temp_row.INDEXER_PARAMETERS) + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/alignment.all.real.rich.vcf.gz" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]) ## take unique vcf file paths
    with open(temp_vcf_path_collection_filename, 'w') as temp_filehandle:
        temp_filehandle.writelines([temp_vcf_path + "\n" for temp_vcf_path in temp_vcf_path_collection])
    return temp_vcf_path_collection_filename



def __write_temp_vcf_dataset_and_sample_and_path_collection_to_file_and_return_filename(dataset_collection_name, temp_vcf_dataset_and_sample_and_path_collection_filename, wildcards, prefix_for_dataset_name_directory="external/DATASET_RNA_EDITING_NAME_DIRECTORY/", prefix_for_dataset_collection_name_directory="external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/"):
    temp_directory = os.path.dirname(temp_vcf_dataset_and_sample_and_path_collection_filename)
    os.makedirs(temp_directory, exist_ok=True)
    temp_vcf_dataset_and_sample_and_path_collection = set([temp_row.DATASET_NAME + "," + temp_row.SAMPLE_NAME + "," + "result/S15_1__get_sample_RNA_editing_sites_v3/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/"  + str(temp_row.INDEXER_PARAMETERS) + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/alignment.all.real.rich.vcf.gz" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]) ## take unique vcf file paths
    with open(temp_vcf_dataset_and_sample_and_path_collection_filename, 'w') as temp_filehandle:
        print('writing ...')
        temp_filehandle.writelines([temp_vcf_dataset_and_sample_and_path + "\n" for temp_vcf_dataset_and_sample_and_path in temp_vcf_dataset_and_sample_and_path_collection])
    return temp_vcf_dataset_and_sample_and_path_collection_filename



def __get_parameters_collection_for_expression_calling(dataset_collection_name, prefix_for_dataset_name_directory="external/DATASET_EXPRESSION_NAME_DIRECTORY/", prefix_for_dataset_collection_name_directory="external/DATASET_EXPRESSION_COLLECTION_NAME_DIRECTORY/"):
    dataset_name_collection = __get_dataset_name_collection(dataset_collection_name, prefix_for_dataset_collection_name_directory=prefix_for_dataset_collection_name_directory)
    final_collection = []
    for temp_dataset_name in dataset_name_collection:
        temp_dataset_DataFrame = pandas.read_csv( prefix_for_dataset_name_directory + temp_dataset_name, sep=",")
        temp_row_collection = [ temp_row for temp_rowindex, temp_row in temp_dataset_DataFrame.iterrows()]
        final_collection = final_collection + temp_row_collection
    print(final_collection)
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

def __get_raw_fastq_path_by_TYPE_and_DATASET_and_SAMPLE_and_RUN(wildcards):
    parent_path = "external/RNA-Seq-with-Run/{DATASET}/{SAMPLE}/{RUN}/RNA/"
    type_end = re.sub("^([^-]+)-.*", "\\1", wildcards.TYPE)
    if type_end == "single":
        return parent_path + "r.fastq.gz"
    elif type_end == "paired":
        return [parent_path + "r1.fastq.gz", parent_path + "r2.fastq.gz"]

## rule S01_1__check_RNA_Seq____step01__check_RNA_Seq_by_fastqc:

rule S01_2__trim_RNA_Seq____step01__trim_RNA_Seq_by_fastp:
    input:
        input_read_filenames=__get_raw_fastq_path_by_TYPE_and_DATASET_and_SAMPLE_and_RUN
    output:
        flag_S01_2__step01=touch("result/S01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/finished.step01__trim_RNA_Seq_by_fastp")
    params:
        raw_fastq_parent_path="external/RNA-Seq-with-Run/{DATASET}/{SAMPLE}/{RUN}/RNA/",
        trim_result_directory="result/S01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/"
    threads:
        config['threads_trimming']
    shell:
        """
        type_end=`echo {wildcards.TYPE} | cut -f 1 -d '-'`
        if [[ $type_end == 'paired' ]]; then
            if [[ "{wildcards.FASTP_PARAMETERS}" == 'none' ]]; then
                ln -s "`readlink -f {params.raw_fastq_parent_path}/r1.fastq.gz`" "{params.trim_result_directory}/r1.fastp.fastq.gz"
                ln -s "`readlink -f {params.raw_fastq_parent_path}/r2.fastq.gz`" "{params.trim_result_directory}/r2.fastp.fastq.gz"
            elif [[ "{wildcards.FASTP_PARAMETERS}" == 'cut-5prime-6-bases-and-base-quality-no-smaller-than-25' ]]; then
                fastp --disable_adapter_trimming --trim_front1 6 --trim_front2 6 --disable_trim_poly_g --disable_quality_filtering --average_qual 25 --disable_length_filtering  -i {params.raw_fastq_parent_path}/r1.fastq.gz -I {params.raw_fastq_parent_path}/r2.fastq.gz -o "{params.trim_result_directory}/r1.fastp.fastq.gz" -O "{params.trim_result_directory}/r2.fastp.fastq.gz" -w {threads} {wildcards.FASTP_PARAMETERS}
            elif [[ "{wildcards.FASTP_PARAMETERS}" == 'base-quality-no-smaller-than-25' ]]; then
                fastp --disable_adapter_trimming --disable_trim_poly_g --disable_quality_filtering --average_qual 25 --disable_length_filtering  -i {params.raw_fastq_parent_path}/r1.fastq.gz -I {params.raw_fastq_parent_path}/r2.fastq.gz -o "{params.trim_result_directory}/r1.fastp.fastq.gz" -O "{params.trim_result_directory}/r2.fastp.fastq.gz" -w {threads} {wildcards.FASTP_PARAMETERS}
            else
                fastp -i {params.raw_fastq_parent_path}/r1.fastq.gz -I {params.raw_fastq_parent_path}/r2.fastq.gz -o "{params.trim_result_directory}/r1.fastp.fastq.gz" -O "{params.trim_result_directory}/r2.fastp.fastq.gz" -w {threads} {wildcards.FASTP_PARAMETERS}
            fi
        elif [[ $type_end == 'single' ]]; then
            if [[ "{wildcards.FASTP_PARAMETERS}" == 'none' ]]; then
                ln -s "`readlink -f {params.raw_fastq_parent_path}/r.fastq.gz`" "{params.trim_result_directory}/r.fastp.fastq.gz"
            elif [[ "{wildcards.FASTP_PARAMETERS}" == 'cut-5prime-6-bases-and-base-quality-no-smaller-than-25' ]]; then
                fastp --disable_adapter_trimming --trim_front1 6 --trim_front2 6 --disable_trim_poly_g --disable_quality_filtering --average_qual 25 --disable_length_filtering  -i {params.raw_fastq_parent_path}/r.fastq.gz -o "{params.trim_result_directory}/r.fastp.fastq.gz" -w {threads} {wildcards.FASTP_PARAMETERS}
            elif [[ "{wildcards.FASTP_PARAMETERS}" == 'base-quality-no-smaller-than-25' ]]; then
                fastp --disable_adapter_trimming --disable_trim_poly_g --disable_quality_filtering --average_qual 25 --disable_length_filtering  -i {params.raw_fastq_parent_path}/r.fastq.gz -o "{params.trim_result_directory}/r.fastp.fastq.gz" -w {threads} {wildcards.FASTP_PARAMETERS}
            else
                fastp -i {params.raw_fastq_parent_path}/r.fastq.gz -o "{params.trim_result_directory}/r.fastp.fastq.gz" -w {threads} {wildcards.FASTP_PARAMETERS}
            fi
        else
            echo "Unsupported type: {wildcards.TYPE}" && exit 1
        fi

        """


## in these series, "TYPE" is of the form "single-50" or "paired-75-60". The numbers are the avgSpotLen (read length) reported by NCBI or some other tools. Currently we ask the users to supply it instead of infer it from fastq files because it is cumbersome in the engineering sense. However, only the first part (single or paired) is used currently. The read length information is just for log purpose only.

## when "{wildcards.INDEXER}" == 'bwa-index-10.1038_nmeth.2330', the "{wildcards.INDEXER_PARAMETER}" is the read length
rule s05_1__index_contig_with_annotation:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}"
    output:
        flag_s05_1=touch("result/s05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/finished")
    params:
        reference_GTF_filename="external/reference.gene.annotation/GENCODE.annotation/{GENCODE_VERSION}/gencode.annotation.gtf",
        reference_UCSC_Table_Browser_knownGene_filename="./external/UCSC.Table.Browser.knownGene.GENCODE/{GENCODE_VERSION}/knownGene",
        contigs_index_directory="result/s05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/",
        contigs_index_prefix="result/s05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/prefix"
    threads:
        config['threads_indexing']
    shell:
        """
        rm -fr "{params.contigs_index_directory}"/*
        mkdir -p "{params.contigs_index_directory}"
        if [[ "{wildcards.INDEXER}" == 'bwa index' ]]; then
            bwa index {wildcards.INDEXER_PARAMETERS} -p "{params.contigs_index_prefix}" "{input.contigs_fasta_filename}"
        elif [[ "{wildcards.INDEXER}" == 'bwa-index-10.1038_nmeth.2330' ]]; then
            cat "{params.reference_UCSC_Table_Browser_knownGene_filename}" | awk '{{print "bin\t"$0}}'  > "{params.contigs_index_directory}"/knownGene_with_dummy_bin_prepended
            perl ./scripts/s05_1__index_contig_with_annotation/getSplicingJunctions_modified.pl "{input.contigs_fasta_filename}" "{params.contigs_index_directory}"/knownGene_with_dummy_bin_prepended  {wildcards.INDEXER_PARAMETERS}  "{params.contigs_index_directory}"/splice_junctions.txt && cat "{input.contigs_fasta_filename}" "{params.contigs_index_directory}"/splice_junctions.txt > "{params.contigs_index_directory}"/combined.fasta
            bwa index -p "{params.contigs_index_prefix}" "{params.contigs_index_directory}"/combined.fasta
        elif [[ "{wildcards.INDEXER}" == 'STAR-expression' ]]; then
            rm -fr {params.contigs_index_prefix}/temp
            rm -fr {params.contigs_index_prefix}
            mkdir -p {params.contigs_index_prefix}
            STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.contigs_index_prefix} --genomeFastaFiles {input.contigs_fasta_filename} --outTmpDir {params.contigs_index_prefix}/temp --sjdbGTFfile {params.reference_GTF_filename} --sjdbOverhang {wildcards.INDEXER_PARAMETERS}
        else
            echo "Indexer {wildcards.INDEXER} not supported." && exit 1;
        fi
        """


rule S15_1__get_sample_RNA_editing_sites_v3__step01__align_sample:
    input:
        flag_S01_2__step01="result/S01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/finished.step01__trim_RNA_Seq_by_fastp",
        flag_s05_1="result/s05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/finished"
    output:
        flag_S15_1__step01=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01__align_sample")
    params:
        contigs_index_prefix="result/s05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/prefix",
        trim_result_directory="result/S01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/",
        alignment_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/"
    threads:
        config['threads_aligning']
    shell:
        """
        type_end=`echo {wildcards.TYPE} | cut -f 1 -d '-'`
        rm -fr "{params.alignment_result_directory}/"*
        mkdir -p "{params.alignment_result_directory}/"
        if [[ "{wildcards.ALIGNER}" == 'bwa-aln-samsepe' ]]; then
            aligner_parameters="{wildcards.ALIGNER_PARAMETERS}"
            bwa_aln_parameters="-h"
            bwa_samse_parameters="-h"
            bwa_sampe_parameters="-h"
            if [[ $aligner_parameters == 'none' ]]; then
                bwa_aln_parameters=""
                bwa_samse_parameters=""
                bwa_sampe_parameters=""
            fi

            if [[ $type_end == single ]]; then
                bwa aln $bwa_aln_parameters -t {threads} -f "{params.alignment_result_directory}/"alignment.sai "{params.contigs_index_prefix}" "{params.trim_result_directory}/r.fastp.fastq.gz"
                bwa samse $bwa_samse_parameters "{params.contigs_index_prefix}" "{params.alignment_result_directory}/"alignment.sai "{params.trim_result_directory}/r.fastp.fastq.gz" | samtools sort -@{threads} -T "{params.alignment_result_directory}/"alignment.bam.temp | samtools view -b -o "{params.alignment_result_directory}/"alignment.bam
            elif [[ $type_end == paired ]]; then
                bwa aln $bwa_aln_parameters -t {threads} -f "{params.alignment_result_directory}/"alignment.1.sai "{params.contigs_index_prefix}" "{params.trim_result_directory}/r1.fastp.fastq.gz"
                bwa aln $bwa_aln_parameters -t {threads} -f "{params.alignment_result_directory}/"alignment.2.sai "{params.contigs_index_prefix}" "{params.trim_result_directory}/r2.fastp.fastq.gz"
                bwa sampe $bwa_sampe_parameters "{params.contigs_index_prefix}" "{params.alignment_result_directory}/"alignment.1.sai "{params.alignment_result_directory}/"alignment.2.sai "{params.trim_result_directory}/r1.fastp.fastq.gz" "{params.trim_result_directory}/r2.fastp.fastq.gz" | samtools sort -@{threads} -T "{params.alignment_result_directory}/"alignment.bam.temp | samtools view -b -o "{params.alignment_result_directory}/"alignment.bam
            else
                echo 'Unsupported read type {wildcards.TYPE}' && exit 1
            fi
        else
            echo 'Unsupported ALIGNER {wildcards.ALIGNER}' && exit 1
        fi
        """


rule s05_2__reformat_dbSNP_vcf:
    input:
        dbSNP_vcf_gz_filename="external/dbSNP.vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/dbSNP.vcf.gz"
    output:
        flag_s05_2=touch("result/s05_2__reformat_dbSNP_vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished"),
        dbSNP_vcf_reformatted_gz_filename="result/s05_2__reformat_dbSNP_vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/dbSNP.vcf.reformatted.gz"
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        zcat {input.dbSNP_vcf_gz_filename} | awk '{{if($1 ~ /^#/){{print $0}} else {{print "chr"$0}} }}' | bgzip -@{threads}  > {output.dbSNP_vcf_reformatted_gz_filename}
        tabix {output.dbSNP_vcf_reformatted_gz_filename}
        """



rule s05_3__index_contig_with_samtools_faidx:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}"
    output:
        flag_s05_3=touch("result/s05_3__index_contig_with_samtools_faidx/{CONTIGS_FASTA_FILENAME}/finished")
    threads:
        1
    shell:
        """
        samtools faidx {input.contigs_fasta_filename}
        """

rule s05_4__create_sequence_dictionary_for_contig_with_picard:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}"
    output:
        flag_s05_4=touch("result/s05_4__create_sequence_dictionary_for_contig_with_picard/{CONTIGS_FASTA_FILENAME}/finished")
    threads:
        1
    shell:
        """
        output_filename=`echo {input.contigs_fasta_filename}|sed 's@\(fa\|fasta\)@dict@'`
        picard CreateSequenceDictionary R={input.contigs_fasta_filename} O=$output_filename
        """





def __get_RUNs_by_DATASET_and_SAMPLE(temp_dataset_name, temp_sample_name, prefix_for_dataset_name_directory="external/DATASET_NAME_DIRECTORY/"):
    temp_dataset_DataFrame = pandas.read_csv( prefix_for_dataset_name_directory + temp_dataset_name, sep=",")
    ## print(temp_dataset_DataFrame)
    temp_run_name_for_the_sample_collection = temp_dataset_DataFrame.loc[temp_dataset_DataFrame["SAMPLE_NAME"] == temp_sample_name, "RUN_NAME"]
    print(temp_run_name_for_the_sample_collection)
    if (len(temp_run_name_for_the_sample_collection) == 0):
        raise ValueError("DATASET: " + temp_dataset_name + " - SAMPLE: " + temp_sample_name + " has no runs available" )
    temp_run_names_list = [ temp_run_name for temp_run_name in temp_run_name_for_the_sample_collection]
    return temp_run_names_list


rule S15_1__get_sample_RNA_editing_sites_v3__step01_2__merge_bams_of_runs_for_the_same_sample:
    input:
        flags_collection_S15_1__step01 = lambda wildcards: ["result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/" + temp_run + "/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01__align_sample" for temp_run in __get_RUNs_by_DATASET_and_SAMPLE(wildcards.DATASET, wildcards.SAMPLE)]
    output:
        flag_S15_1__step01_2=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01_2__merge_bams_of_runs_for_the_same_sample"),
        alignment_merged_filename="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/alignment.merged.bam"
    params:
        alignment_bams_prefix = "result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/",
        alignment_bams_suffix = "/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/alignment.bam",
        alignment_bams_runs = lambda wildcards: __get_RUNs_by_DATASET_and_SAMPLE(wildcards.DATASET, wildcards.SAMPLE),
        number_of_alignments_to_merge = lambda wildcards: len(__get_RUNs_by_DATASET_and_SAMPLE(wildcards.DATASET, wildcards.SAMPLE)),
        merge_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/",
        picard_Xmx=config["picard_Xmx"]
    threads:
        config['threads_merging_bams']
    shell:
        """
        if [[ {params.number_of_alignments_to_merge} > 1 ]]; then
            final_INPUT_literal=''
            for item in {params.alignment_bams_runs:q}
            do
                final_INPUT_literal=$final_INPUT_literal"INPUT="{params.alignment_bams_prefix}$item{params.alignment_bams_suffix}" "
            done
            picard MergeSamFiles {params.picard_Xmx} $final_INPUT_literal OUTPUT={output.alignment_merged_filename} SO=coordinate VALIDATION_STRINGENCY=SILENT
        else ## == 1
            ln -s `readlink -f {params.alignment_bams_prefix}{params.alignment_bams_runs}{params.alignment_bams_suffix}` {output.alignment_merged_filename}
        fi
        """



rule S15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part01__sort_bam:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}",
        flag_s05_3="result/s05_3__index_contig_with_samtools_faidx/{CONTIGS_FASTA_FILENAME}/finished",
        flag_s05_4="result/s05_4__create_sequence_dictionary_for_contig_with_picard/{CONTIGS_FASTA_FILENAME}/finished",
        flag_S15_1__step01_2="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01_2__merge_bams_of_runs_for_the_same_sample"
    output:
        flag_S15_1__step02__part01=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part01__sort_bam")
    params:
        alignment_merged_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/",
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        picard_Xmx=config["picard_Xmx"]
    threads:
        config['threads_calling_variants']
    shell:
        """
        mkdir -p "{params.calling_result_directory}/"
        picard SortSam {params.picard_Xmx} INPUT="{params.alignment_merged_directory}/alignment.merged.bam" OUTPUT="{params.calling_result_directory}/alignment.sorted.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
        """

rule S15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part02__label_bam:
    input:
        flag_S15_1__step02__part01="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part01__sort_bam"
    output:
        flag_S15_1__step02__part02=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part02__label_bam")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        picard_Xmx=config["picard_Xmx"]
    threads:
        config['threads_calling_variants']
    shell:
        """
        picard AddOrReplaceReadGroups  {params.picard_Xmx}  I="{params.calling_result_directory}/alignment.sorted.bam" O="{params.calling_result_directory}/alignment.sorted.withRG.bam" SO=coordinate RGID="{wildcards.SAMPLE}" RGLB="{wildcards.SAMPLE}" RGPL=illumina RGPU=illumina RGSM="{wildcards.SAMPLE}" VALIDATION_STRINGENCY=SILENT
        """


rule S15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part03__mark_duplicates:
    input:
        flag_S15_1__step02__part02="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part02__label_bam"
    output:
        flag_S15_1__step02__part03=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part03__mark_duplicates")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        picard_Xmx=config["picard_Xmx"]
    threads:
        config['threads_calling_variants']
    shell:
        """
        picard MarkDuplicates  {params.picard_Xmx}  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=2000 INPUT="{params.calling_result_directory}/alignment.sorted.withRG.bam" OUTPUT="{params.calling_result_directory}/alignment.sorted.withRG.dedup.bam" METRICS_FILE="{params.calling_result_directory}/alignment.sorted.withRG.dedup.matrix" CREATE_INDEX=true, VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
        """


rule S15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part04__correct_coordinates:
    input:
        flag_S15_1__step02__part03="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part03__mark_duplicates"
    output:
        flag_S15_1__step02__part04=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part04__correct_coordinates")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        convertCoordinates_Xmx=config["convertCoordinates_Xmx"],
        picard_Xmx=config["picard_Xmx"]
    threads:
        config['threads_calling_variants']
    shell:
        """
        samtools view -h "{params.calling_result_directory}/alignment.sorted.withRG.dedup.bam" > "{params.calling_result_directory}/alignment.sorted.withRG.dedup.sam"
        java {params.convertCoordinates_Xmx} -classpath tools/convertCoordinates_classpath convertCoordinates < "{params.calling_result_directory}/alignment.sorted.withRG.dedup.sam" | samtools view -b -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bam"
        rm -fr "{params.calling_result_directory}/alignment.sorted.withRG.dedup.sam"
        samtools view -Sq 20 -h "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bam" -b -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.bam"
        picard SortSam  {params.picard_Xmx}  INPUT="{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.bam" OUTPUT="{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=YES
        samtools view -h "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.bam" | grep -v -Pe "SN:chr[1234567890XY]+-" | samtools view -@1 -b -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
        samtools index -@1 "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
        """



rule S15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part05__locally_realign_reads:
    input:
        flag_S15_1__step02__part04="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part04__correct_coordinates",
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}",
        flag_s05_2="result/s05_2__reformat_dbSNP_vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished",
        dbSNP_vcf_reformatted_gz_filename="result/s05_2__reformat_dbSNP_vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/dbSNP.vcf.reformatted.gz"
    output:
        flag_S15_1__step02__part05=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part05__locally_realign_reads")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        GATK_RealignerTargetCreator_Xmx=config["GATK_RealignerTargetCreator_Xmx"],
        GATK_IndelRealigner_Xmx=config["GATK_IndelRealigner_Xmx"],
        GATK_BaseRecalibrator_Xmx=config["GATK_BaseRecalibrator_Xmx"],
        GATK_PrintReads_Xmx=config["GATK_PrintReads_Xmx"]
    threads:
        config['threads_calling_variants']
    shell:
        """
        java {params.GATK_RealignerTargetCreator_Xmx} -jar ./tools/{wildcards.VARIANT_CALLER}/GenomeAnalysisTK.jar -nt 1 -R "{input.contigs_fasta_filename}" -T RealignerTargetCreator -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.realn.intervals" -I "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam" -U ALLOW_N_CIGAR_READS
        java {params.GATK_IndelRealigner_Xmx} -jar ./tools/{wildcards.VARIANT_CALLER}/GenomeAnalysisTK.jar -R "{input.contigs_fasta_filename}" -T IndelRealigner -targetIntervals "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.realn.intervals" -I "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam" -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.realn.bam" -U ALLOW_N_CIGAR_READS
        java {params.GATK_BaseRecalibrator_Xmx} -jar ./tools/{wildcards.VARIANT_CALLER}/GenomeAnalysisTK.jar -R "{input.contigs_fasta_filename}" -T BaseRecalibrator -I "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.realn.bam" -knownSites {input.dbSNP_vcf_reformatted_gz_filename} -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal_data.grp" -U ALLOW_N_CIGAR_READS
        java {params.GATK_PrintReads_Xmx} -jar ./tools/{wildcards.VARIANT_CALLER}/GenomeAnalysisTK.jar -R "{input.contigs_fasta_filename}" -T PrintReads -I "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.realn.bam" -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam" -BQSR "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal_data.grp" -U ALLOW_N_CIGAR_READS
        """

rule S15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part06__really_call_variants:
    input:
        flag_S15_1__step02__part05="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part05__locally_realign_reads",
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}",
        flag_s05_2="result/s05_2__reformat_dbSNP_vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished",
        dbSNP_vcf_reformatted_gz_filename="result/s05_2__reformat_dbSNP_vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/dbSNP.vcf.reformatted.gz"
    output:
        flag_S15_1__step02__part06=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part06__really_call_variants"),
        alignment_bcf_filename="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/alignment.bcf"
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        GATK_UnifiedGenotyper_Xmx=config["GATK_UnifiedGenotyper_Xmx"],
        GATK_VariantFiltration_Xmx=config["GATK_VariantFiltration_Xmx"]
    threads:
        config['threads_calling_variants']
    shell:
        """
        java {params.GATK_UnifiedGenotyper_Xmx} -jar ./tools/{wildcards.VARIANT_CALLER}/GenomeAnalysisTK.jar -nt 1 -R "{input.contigs_fasta_filename}" -T UnifiedGenotyper -I "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam" -o  "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.gatk.raw.vcf.gz" --dbsnp {input.dbSNP_vcf_reformatted_gz_filename} -stand_call_conf 0 -stand_emit_conf 0 -U ALLOW_N_CIGAR_READS
        java {params.GATK_VariantFiltration_Xmx} -jar ./tools/{wildcards.VARIANT_CALLER}/GenomeAnalysisTK.jar -nt 1 -R "{input.contigs_fasta_filename}" -T VariantFiltration  --variant "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.gatk.raw.vcf.gz" -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.gatk.raw.dummy.filtered.vcf.gz"
        bcftools view -o "{params.calling_result_directory}/alignment.bcf" -Ob "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.gatk.raw.dummy.filtered.vcf.gz"
        """


rule S15_1__get_sample_RNA_editing_sites_v3____step07__apply_complex_filter____part01__remove_dbSNP_but_keep_cDNA_flagged_dbSNP:
    input:
        flag_S15_1__step02="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part06__really_call_variants",
        alignment_bcf_filename="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/alignment.bcf"
    output:
        flag_S15_1__step07__part01=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part01__remove_dbSNP_but_keep_cDNA_flagged_dbSNP")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        complex_filter_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/",
        flagged_cDNA_only_dbSNP_bed_filename="external/UCSC.Table.Browser.dbSNP/{dbSNP_VERSION}/flagged.cDNA.only/dbSNP.bed"
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        rm -fr "{params.complex_filter_directory}"/*
        mkdir -p "{params.complex_filter_directory}"
        perl ./scripts/S15_1__get_sample_RNA_editing_sites_v3/step07__apply_complex_filter/complex_filter_1/Convert_VCF.pl <(bcftools view -Ov "{params.calling_result_directory}/alignment.bcf") "{params.complex_filter_directory}"/alignment.con.vcf
        perl ./scripts/S15_1__get_sample_RNA_editing_sites_v3/step07__apply_complex_filter/complex_filter_1/ref_filter.pl "{params.flagged_cDNA_only_dbSNP_bed_filename}" "{params.complex_filter_directory}"/alignment.con.vcf  <(bcftools view -Ov "{params.calling_result_directory}/alignment.bcf") "{params.complex_filter_directory}"/alignment.ref.vcf
        """

rule S15_1__get_sample_RNA_editing_sites_v3____step07__apply_complex_filter____part02__remove_mismatch_at_first6bp:
    input:
        flag_S15_1__step07__part01="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part01__remove_dbSNP_but_keep_cDNA_flagged_dbSNP"
    output:
        flag_S15_1__step07__part02=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part02__remove_mismatch_at_first6bp")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        complex_filter_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/"
    threads:
        config['threads_auxiliary_processing_parallel']
    shell:
        """
        samtools index -@{threads} "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.realn.bam"
        perl ./scripts/S15_1__get_sample_RNA_editing_sites_v3/step07__apply_complex_filter/complex_filter_1/threads_rmMismatch.pl "{params.complex_filter_directory}"/alignment.ref.vcf "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.realn.bam" "{params.complex_filter_directory}"/alignment.rem.vcf {threads}
        """


rule S15_1__get_sample_RNA_editing_sites_v3____step07__apply_complex_filter____part03__separate_Alu_and_others:
    input:
        flag_S15_1__step07__part02="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part02__remove_mismatch_at_first6bp"
    output:
        flag_S15_1__step07__part03=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part03__separate_Alu_and_others")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        complex_filter_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/",
        repeatmasker_repFamily_Alu_bed_filename="external/UCSC.Table.Browser.repeatmasker/repFamily.Alu/repeatmasker.bed"
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        perl ./scripts/S15_1__get_sample_RNA_editing_sites_v3/step07__apply_complex_filter/complex_filter_1/Alu_filter.pl "{params.repeatmasker_repFamily_Alu_bed_filename}" "{params.complex_filter_directory}"/alignment.rem.vcf "{params.complex_filter_directory}"/alignment.Alu.vcf "{params.complex_filter_directory}"/alignment.others.vcf
        """



rule S15_1__get_sample_RNA_editing_sites_v3____step07__apply_complex_filter____part04__filter_against_sites_in_simple_repeats:
    input:
        flag_S15_1__step07__part03="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part03__separate_Alu_and_others"
    output:
        flag_S15_1__step07__part04=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part04__filter_against_sites_in_simple_repeats")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        complex_filter_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/",
        repeatmasker_repFamily_Simple_repeat_bed_filename="external/UCSC.Table.Browser.repeatmasker/repFamily.Simple_repeat/repeatmasker.bed"
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        perl ./scripts/S15_1__get_sample_RNA_editing_sites_v3/step07__apply_complex_filter/complex_filter_1/FreqSimple.pl "{params.repeatmasker_repFamily_Simple_repeat_bed_filename}" "{params.complex_filter_directory}"/alignment.others.vcf  "{params.complex_filter_directory}"/alignment.others.sim.vcf
        """



rule S15_1__get_sample_RNA_editing_sites_v3____step07__apply_complex_filter____part05__filter_against_sites_near_splice_junctions_or_homopolymers:
    input:
        flag_S15_1__step07__part04="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part04__filter_against_sites_in_simple_repeats"
    output:
        flag_S15_1__step07__part05=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part05__filter_against_sites_near_splice_junctions_or_homopolymers")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        complex_filter_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/",
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}",
        contigs_index_directory="result/s05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/"
    threads:
        config['threads_auxiliary_processing_parallel']
    shell:
        """
        perl ./scripts/S15_1__get_sample_RNA_editing_sites_v3/step07__apply_complex_filter/complex_filter_1/threads_rmSJandHomo.pl "{params.contigs_fasta_filename}" "{params.complex_filter_directory}"/alignment.others "{params.contigs_index_directory}"/knownGene_with_dummy_bin_prepended {threads}
        """



rule S15_1__get_sample_RNA_editing_sites_v3____step07__apply_complex_filter____part06__filter_against_multimapped_reads_by_blat:
    input:
        flag_S15_1__step07__part05="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part05__filter_against_sites_near_splice_junctions_or_homopolymers"
    output:
        flag_S15_1__step07__part06=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part06__filter_against_multimapped_reads_by_blat")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        complex_filter_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/",
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}",
        contigs_index_directory="result/s05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/"
    threads:
        config['threads_auxiliary_processing_parallel']
    shell:
        """
        perl ./scripts/S15_1__get_sample_RNA_editing_sites_v3/step07__apply_complex_filter/complex_filter_1/threads_BlatCandidates.pl "{params.contigs_fasta_filename}" "{params.complex_filter_directory}"/alignment.others.rmSJandHomo.txt "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.realn.bam" "{params.complex_filter_directory}"/alignment.others.blat.vcf {threads}
        """


rule S15_1__get_sample_RNA_editing_sites_v3____step07__apply_complex_filter____part07__separate_RepNOTAlu_and_nonRep:
    input:
        flag_S15_1__step07__part06="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part06__filter_against_multimapped_reads_by_blat"
    output:
        flag_S15_1__step07__part07=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part07__separate_RepNOTAlu_and_nonRep")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        complex_filter_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/",
        repeatmasker_repFamily_not_Alu_bed_filename="external/UCSC.Table.Browser.repeatmasker/repFamily.not.Alu/repeatmasker.bed"
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        perl ./scripts/S15_1__get_sample_RNA_editing_sites_v3/step07__apply_complex_filter/complex_filter_1/nonAlu_filter_new.pl "{params.repeatmasker_repFamily_not_Alu_bed_filename}" "{params.complex_filter_directory}"/alignment.others.blat.vcf "{params.complex_filter_directory}"/alignment.RepNOTAlu.vcf "{params.complex_filter_directory}"/alignment.nonRep.vcf
        """


rule S15_1__get_sample_RNA_editing_sites_v3____step07__apply_complex_filter____part08__reformat_data_as_standard_vcf:
    input:
        flag_S15_1__step07__part07="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part07__separate_RepNOTAlu_and_nonRep"
    output:
        flag_S15_1__step07__part08=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part08__reformat_data_as_standard_vcf")
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        complex_filter_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/"
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        for subset in Alu RepNOTAlu nonRep
        do
            ## create VCF that can be accepted by bedtools
            echo "##fileformat=VCF" > "{params.complex_filter_directory}/alignment.$subset".real.vcf
            echo "#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO" | tr ',' '\t'  >> "{params.complex_filter_directory}/alignment.$subset".real.vcf
            cat "{params.complex_filter_directory}/alignment.$subset".vcf | awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,$4,$5,$6,".","."}}' >> "{params.complex_filter_directory}/alignment.$subset".real.vcf
        done
        """

rule B15_1__get_sample_RNA_editing_sites_v3____step07__apply_complex_filter_1____part08__reformat_data_as_standard_vcf:
    input:
        lambda wildcards: ["result/S15_1__get_sample_RNA_editing_sites_v3/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/" + str(temp_row.INDEXER_PARAMETERS) + "/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part08__reformat_data_as_standard_vcf" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_B15_1__step10=touch("result/B15_1__get_sample_RNA_editing_sites_v3/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part08__reformat_data_as_standard_vcf")





rule B15_1__get_sample_RNA_editing_sites_v3____step07__apply_complex_filter_1____part09__reformat_data_as_standard_rich_vcf:
    input:
        lambda wildcards: ["result/S15_1__get_sample_RNA_editing_sites_v3/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/" + str(temp_row.INDEXER_PARAMETERS) + "/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part09__reformat_data_as_standard_rich_vcf" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_B15_1__step07__part09=touch("result/B15_1__get_sample_RNA_editing_sites_v3/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part09__reformat_data_as_standard_rich_vcf")


rule S16_1__concatenate_RNA_editing_site_from_a_dataset_collection:
    input:
        flag_B15_1__step10="result/B15_1__get_sample_RNA_editing_sites_v3/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part09__reformat_data_as_standard_rich_vcf"
    output:
        flag_S16_1=touch("result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished"),
        merged_vcf_gz_filename="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/merged.vcf.gz",
        merged_variant_only_vcf_gz_filename="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/merged.variant.only.vcf.gz"
    params:
        temp_vcf_path_collection = lambda wildcards: __write_temp_vcf_path_collection_to_file_and_return_filename(dataset_collection_name=wildcards.DATASET_COLLECTION_NAME, temp_vcf_path_collection_filename= "result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/vcf_path_to_merge.txt", wildcards=wildcards),
        result_directory="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/"
    threads:
        config['threads_merging_vcfs']
    shell:
        """
        rm -fr {params.result_directory}/*vcf.gz
        mkdir -p {params.result_directory}
        bcftools merge --threads {threads} -o {params.result_directory}/merged.vcf.gz -Oz --file-list {params.result_directory}/vcf_path_to_merge.txt
        tabix {params.result_directory}/merged.vcf.gz
        bcftools view --threads {threads} --samples nonexistent --force-samples -o {params.result_directory}/merged.variant.only.vcf.gz -Oz {params.result_directory}/merged.vcf.gz
        tabix {params.result_directory}/merged.variant.only.vcf.gz
        """





rule x16_2__annotate_merged_sites_from_a_dataset_collection____step01__annotate_with_snpEff____patch01__build_database:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}",
        reference_GTF_filename="external/reference.gene.annotation/GENCODE.annotation/{GENCODE_VERSION}/gencode.annotation.gtf"
    output:
        flag_x05_1_step11_patch01=touch("result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.patch01__build_database")
    params:
        snpEff_data_directory="result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/snpEff.data",
        snpEff_config_filename="result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/snpEff.config",
        snpEff_database_name="{CONTIGS_FASTA_FILENAME}.GENCODE.{GENCODE_VERSION}",
        snpEff_Xmx=config['snpEff_Xmx']
    threads:
        config['threads_indexing']
    shell:
        """
        rm -fr {params.snpEff_data_directory}/*
        mkdir -p {params.snpEff_data_directory}/{params.snpEff_database_name}
        echo data.dir = `readlink -f {params.snpEff_data_directory}` > {params.snpEff_config_filename}
        echo {wildcards.CONTIGS_FASTA_FILENAME} with GENCODE {wildcards.GENCODE_VERSION},Version 1 >> {params.snpEff_config_filename}
        echo {params.snpEff_database_name}.genome : {params.snpEff_database_name} >> {params.snpEff_config_filename}
        ln -s `readlink -f {input.reference_GTF_filename}` {params.snpEff_data_directory}/{params.snpEff_database_name}/genes.gtf
        ln -s `readlink -f {input.contigs_fasta_filename}` {params.snpEff_data_directory}/{params.snpEff_database_name}/sequences.fa
        snpEff build {params.snpEff_Xmx} -gtf22 -config {params.snpEff_config_filename} -verbose {params.snpEff_database_name}
        """


rule S16_2__annotate_merged_sites_from_a_dataset_collection____step01__annotate_with_snpEff:
    input:
        flag_S16_1="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished",
        merged_variant_only_vcf_gz_filename="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/merged.variant.only.vcf.gz",
        flag_x05_1_step11_patch01="result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.patch01__build_database"
    output:
        flag_S16_2_step01=touch("result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/finished.step01__annotate_with_snpEff"),
        merged_variant_only_snpEff_vcf_gz_filename="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/merged.variant.only.snpEff.vcf.gz"
    params:
        snpEff_data_directory="result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/snpEff.data",
        snpEff_config_filename="result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/snpEff.config",
        snpEff_database_name="{CONTIGS_FASTA_FILENAME}.GENCODE.{GENCODE_VERSION}",
        snpEff_Xmx=config['snpEff_Xmx'],
        annotation_directory="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/"
    threads:
        config["threads_annotating"]
    shell:
        """
        snpEff ann {params.snpEff_Xmx} -lof -config {params.snpEff_config_filename} {params.snpEff_database_name} {input.merged_variant_only_vcf_gz_filename} | bcftools view -o {output.merged_variant_only_snpEff_vcf_gz_filename} -Oz
        tabix {output.merged_variant_only_snpEff_vcf_gz_filename}
        """



rule S16_2__annotate_merged_sites_from_a_dataset_collection____step02__bcftools_isec_with_outer_vcf:
    input:
        flag_S16_1="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished",
        merged_variant_only_vcf_gz_filename="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/merged.variant.only.vcf.gz",
        outer_vcf_filename="external/outer_vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/outer.VCF",
        outer_vcf_tbi_filename="external/outer_vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/outer.VCF.tbi",
        flag_x05_1_step11_patch01="result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.patch01__build_database"
    output:
        flag_S16_2_step02=touch("result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/finished.step02__bcftools_isec_with_outer_vcf"),
        bcftools_isec_outer_vcf_result_vcf_filename="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/bcftools.isec.outer.vcf.result.vcf.gz"
    params:
        result_directory="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/"
    threads:
        config["threads_annotating"]
    shell:
        """
        rm -fr {params.result_directory}/*
        mkdir -p {params.result_directory}
        bcftools isec --threads {threads} {input.merged_variant_only_vcf_gz_filename} {input.outer_vcf_filename} -n =2 -w 2 -o {output.bcftools_isec_outer_vcf_result_vcf_filename} -Oz
        tabix {output.bcftools_isec_outer_vcf_result_vcf_filename}
        """


rule S16_2__annotate_merged_sites_from_a_dataset_collection____step03__compute_tables_and_stats_based_on_snpEff_annotation:
    input:
        flag_S16_1="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished",
        merged_vcf_gz_filename="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/merged.vcf.gz",
        flag_S16_2_step01="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/finished.step01__annotate_with_snpEff",
        merged_variant_only_snpEff_vcf_gz_filename="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/merged.variant.only.snpEff.vcf.gz"
    output:
        flag_S16_2_step03=touch("result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/finished.step03__compute_tables_and_stats_based_on_snpEff_annotation"),
        merged_variant_only_snpEff_reformatted_dt_filename="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/merged.variant.only.snpEff.reformatted.dt.txt",
        merged_variant_only_snpEff_ANN_split_dt_filename="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/merged.variant.only.snpEff.ANN.split.dt.txt",
        merged_variant_only_snpEff_ANN_single_match_dt_filename="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/merged.variant.only.snpEff.ANN.single.match.dt.txt",
        merged_variant_only_snpEff_event_summary_dt_filename="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/merged.variant.only.snpEff.event.summary.dt.txt",
        merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_filename="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/merged.vcf.reformatted.with.snpEff.ANN.split.annotation.dt.txt",
        merged_vcf_reformatted_single_event_only_melt_dt_filename="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/merged.vcf.reformatted.single.event.only.melt.dt.txt",
        merged_vcf_reformatted_single_event_only_melt_with_snpEff_ANN_split_annotation_dt_filename="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/merged.vcf.reformatted.single.event.only.melt.with.snpEff.ANN.split.annotation.dt.txt",
        merged_vcf_reformatted_AEI_per_sample_dt_filename="result/S16_2__annotate_merged_sites_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/AEI.per.sample.dt.txt"
    threads:
        config["threads_auxiliary_processing"]
    script:
        "./scripts/S16_2__annotate_merged_sites_from_a_dataset_collection/step03__compute_tables_and_stats_based_on_snpEff_annotation/run.R"



rule S16_3__get_RNA_editing_site_long_table_from_a_dataset_collection:
    input:
        flag_B15_1__step10="result/B15_1__get_sample_RNA_editing_sites_v3/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part09__reformat_data_as_standard_rich_vcf"
    output:
        flag_S16_3=touch("result/S16_3__get_RNA_editing_site_long_table_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished"),
        merged_long_table_dt_txt_gz_filename="result/S16_3__get_RNA_editing_site_long_table_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/merged.long.table.dt.txt.gz"
    params:
        temp_vcf_dataset_and_sample_and_path_collection = lambda wildcards: __write_temp_vcf_dataset_and_sample_and_path_collection_to_file_and_return_filename(dataset_collection_name=wildcards.DATASET_COLLECTION_NAME, temp_vcf_dataset_and_sample_and_path_collection_filename= "result/S16_3__get_RNA_editing_site_long_table_from_a_dataset_collection/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/vcf_dataset_and_sample_and_path_to_merge.txt", wildcards=wildcards),
        temp_result_directory="result/S16_3__get_RNA_editing_site_long_table_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/"
    threads:
        config['threads_merging_vcfs']
    shell:
        """
        rm -fr {params.temp_result_directory}/temp.result.txt && mkdir -p {params.temp_result_directory}
        echo ID,SUBSET,SAMPLE,AC,AN,AF | tr ',' '\t' > {params.temp_result_directory}/temp.result.txt
        cat {params.temp_vcf_dataset_and_sample_and_path_collection} | while read line; do
            echo `date` Processing $line
            dataset=`echo $line | cut -f 1 -d ,`
            sample=`echo $line | cut -f 2 -d ,`
            vcf=`echo $line | cut -f 3 -d ,`
            bcftools query -f "%ID\t%INFO/SUBSET\t$sample\t[%AC\t%AN\t%AF]\n" $vcf >> {params.temp_result_directory}/temp.result.txt
        done
        pigz --processes {threads} --stdout {params.temp_result_directory}/temp.result.txt > {output.merged_long_table_dt_txt_gz_filename}
        """


rule S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome:
    input:
        flag_B15_1__step10="result/B15_1__get_sample_RNA_editing_sites_v3/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part09__reformat_data_as_standard_rich_vcf"
    output:
        flag_S17_1=touch("result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/finished"),
        merged_vcf_gz_filename="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/merged.vcf.gz",
        merged_variant_only_vcf_gz_filename="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/merged.variant.only.vcf.gz"
    params:
        vcf_path_to_merge_filename = lambda wildcards: __write_temp_vcf_path_collection_to_file_and_return_filename(dataset_collection_name=wildcards.DATASET_COLLECTION_NAME, temp_vcf_path_collection_filename= "result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + wildcards.REGIONS + "/vcf_path_to_merge.txt", wildcards=wildcards),
        result_directory="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}"
    threads:
        config['threads_merging_vcfs']
    shell:
        """
        rm -fr {params.result_directory}/*vcf.gz
        mkdir -p {params.result_directory}
        bcftools merge --regions {wildcards.REGIONS} --threads {threads} -o {params.result_directory}/merged.vcf.gz -Oz --file-list {params.vcf_path_to_merge_filename}
        tabix {params.result_directory}/merged.vcf.gz
        bcftools view --threads {threads} --samples nonexistent --force-samples -o {params.result_directory}/merged.variant.only.vcf.gz -Oz {params.result_directory}/merged.vcf.gz
        tabix {params.result_directory}/merged.variant.only.vcf.gz
        """


rule S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome____step01__annotate_with_snpEff:
    input:
        flag_S17_1="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/finished",
        merged_variant_only_vcf_gz_filename="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/merged.variant.only.vcf.gz",
        flag_x16_2_step01_patch01="result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.patch01__build_database"
    output:
        flag_S17_2_step01=touch("result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/finished.step01__annotate_with_snpEff"),
        merged_variant_only_snpEff_vcf_gz_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/merged.variant.only.snpEff.vcf.gz"
    params:
        snpEff_data_directory="result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/snpEff.data",
        snpEff_config_filename="result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/snpEff.config",
        snpEff_database_name="{CONTIGS_FASTA_FILENAME}.GENCODE.{GENCODE_VERSION}",
        snpEff_Xmx=config['snpEff_Xmx'],
        annotation_directory="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/"
    threads:
        config["threads_annotating"]
    shell:
        """
        snpEff ann {params.snpEff_Xmx} -lof -config {params.snpEff_config_filename} {params.snpEff_database_name} {input.merged_variant_only_vcf_gz_filename} | bcftools view -o {output.merged_variant_only_snpEff_vcf_gz_filename} -Oz
        tabix {output.merged_variant_only_snpEff_vcf_gz_filename}
        """


rule S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome___step03__compute_tables_and_stats_based_on_snpEff_annotation:
    input:
        flag_S17_1="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/finished",
        merged_vcf_gz_filename="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/merged.vcf.gz",
        flag_S17_2_step01="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/finished.step01__annotate_with_snpEff",
        merged_variant_only_snpEff_vcf_gz_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/merged.variant.only.snpEff.vcf.gz"
    output:
        flag_S17_2_step03=touch("result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/finished.step03__compute_tables_and_stats_based_on_snpEff_annotation"),
        merged_variant_only_snpEff_reformatted_dt_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/merged.variant.only.snpEff.reformatted.dt.txt",
        merged_variant_only_snpEff_ANN_split_dt_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/merged.variant.only.snpEff.ANN.split.dt.txt",
        merged_variant_only_snpEff_ANN_single_match_dt_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/merged.variant.only.snpEff.ANN.single.match.dt.txt",
        merged_variant_only_snpEff_event_summary_dt_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/merged.variant.only.snpEff.event.summary.dt.txt",
        merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/merged.vcf.reformatted.with.snpEff.ANN.split.annotation.dt.txt",
        merged_vcf_reformatted_single_event_only_melt_dt_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/merged.vcf.reformatted.single.event.only.melt.dt.txt",
        merged_vcf_reformatted_single_event_only_melt_with_snpEff_ANN_split_annotation_dt_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/merged.vcf.reformatted.single.event.only.melt.with.snpEff.ANN.split.annotation.dt.txt",
        AEI_per_sample_dt_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/snpEff/AEI.per.sample.dt.txt"
    threads:
        config["threads_auxiliary_processing"]
    script:
        "./scripts/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/step03__compute_tables_and_stats_based_on_snpEff_annotation/run.R"




rule S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome____step04__find_overlaps_with_known_DNA_variants:
    input:
        flag_S17_1="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/finished",
        merged_variant_only_vcf_gz_filename="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/merged.variant.only.vcf.gz",
        outer_vcf_filename="external/outer_vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/outer.VCF",
        outer_vcf_tbi_filename="external/outer_vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/outer.VCF.tbi"
    output:
        flag_S17_2_step04=touch("result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/finished.step04__find_overlaps_with_known_DNA_variants"),
        merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/merged.variant.only.bcftools.isec.outer.vcf.result.vcf.gz"
    params:
        annotation_directory="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/"
    threads:
        config["threads_annotating"]
    shell:
        """
        rm -fr {params.annotation_directory}/*
        mkdir -p {params.annotation_directory}
        bcftools isec --threads {threads} {input.merged_variant_only_vcf_gz_filename} {input.outer_vcf_filename} -n =2 -w 2 -o {output.merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename} -Oz
        tabix {output.merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename}
        """



rule S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome____step05__find_overlaps_with_known_DNA_variants_with_collapse_all:
    input:
        flag_S16_1="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/finished",
        merged_variant_only_vcf_gz_filename="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/merged.variant.only.vcf.gz",
        outer_vcf_filename="external/outer_vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/outer.VCF",
        outer_vcf_tbi_filename="external/outer_vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/outer.VCF.tbi"
    output:
        flag_S17_2_step05=touch("result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/collapse_all/finished.step05__find_overlaps_with_known_DNA_variants_with_collapse_all"),
        merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/collapse_all/merged.variant.only.bcftools.isec.outer.vcf.result.vcf.gz"
    params:
        annotation_directory="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/collapse_all/"
    threads:
        config["threads_annotating"]
    shell:
        """
        rm -fr {params.annotation_directory}/*
        mkdir -p {params.annotation_directory}
        bcftools isec --collapse all --threads {threads} {input.merged_variant_only_vcf_gz_filename} {input.outer_vcf_filename} -n =2 -w 2 -o {output.merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename} -Oz
        tabix {output.merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename}
        """



rule S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome____step06__find_overlaps_with_known_DNA_variants_with_collapse_all_and_keep_self_vcf:
    input:
        flag_S16_1="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/finished",
        merged_variant_only_vcf_gz_filename="result/S17_1__concatenate_RNA_editing_site_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/merged.variant.only.vcf.gz",
        outer_vcf_filename="external/outer_vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/outer.VCF",
        outer_vcf_tbi_filename="external/outer_vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/outer.VCF.tbi"
    output:
        flag_S17_2_step06=touch("result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/collapse_all_and_keep_self_vcf/finished.step06__find_overlaps_with_known_DNA_variants_with_collapse_all_and_keep_self_vcf"),
        merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_self_vcf_gz_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/collapse_all_and_keep_self_vcf/merged.variant.only.bcftools.isec.outer.vcf.result.self.vcf.gz",
        merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_outer_vcf_gz_filename="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/collapse_all_and_keep_self_vcf/merged.variant.only.bcftools.isec.outer.vcf.result.outer.vcf.gz"
    params:
        annotation_directory="result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{REGIONS}/bcftools.isec.with.outer.vcf/{OUTER_VCF_NAME}/{OUTER_VCF_SUBSET}/collapse_all_and_keep_self_vcf/"
    threads:
        config["threads_annotating"]
    shell:
        """
        rm -fr {params.annotation_directory}/*
        mkdir -p {params.annotation_directory}/raw.output/
        bcftools isec --collapse all --threads {threads} {input.merged_variant_only_vcf_gz_filename} {input.outer_vcf_filename} -n=2 -p {params.annotation_directory}/raw.output/ -Oz
        ln -s -r {params.annotation_directory}/raw.output/0000.vcf.gz {output.merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_self_vcf_gz_filename}
        ln -s -r {params.annotation_directory}/raw.output/0001.vcf.gz {output.merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_outer_vcf_gz_filename}
        tabix {output.merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_self_vcf_gz_filename}
        tabix {output.merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_outer_vcf_gz_filename}
        """



def __get_uniform_regions(chromosomes, bin_size):
    ##currently only support chromosomes=='basic'
    assert chromosomes == 'basic' or chromosomes == 'basic_without_Y', "`chromosomes` must be 'basic' or 'basic_without_Y', not " + str(chromosomes)
    bin_size_integer = int(bin_size)
    final_collection = []

    non_numeric_chromosome_indices_collection = ["X", "Y"]
    if chromosomes == 'basic_without_Y':
        non_numeric_chromosome_indices_collection = ["X"]
    for temp_chromosome_index in [str(temp_i) for temp_i in range(1, 23)] + non_numeric_chromosome_indices_collection:
        for temp_bin_start_closed in range(1, 250000000, bin_size_integer):
            temp_bin_end_closed = temp_bin_start_closed + bin_size_integer - 1
            final_collection.append(("chr" + temp_chromosome_index, temp_bin_start_closed, temp_bin_end_closed))
    return final_collection

rule S18_1__combine_annotations____step01__combine_merged_vcf_reformatted_single_event_only_melt_dt_filename:
    input:
        merged_vcf_reformatted_single_event_only_melt_dt_filename_collection=lambda wildcards: ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/snpEff/merged.vcf.reformatted.single.event.only.melt.dt.txt" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    output:
        flag_S18_1_step01=touch("result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/{CHROMOSOMES}/{BIN_SIZE}/finished.step01__combine_merged_vcf_reformatted_single_event_only_melt_dt_filename"),
        combined_merged_vcf_reformatted_single_event_only_melt_dt_filename="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/snpEff/{COMPLEX_FILTER_PARAMETERS}/{CHROMOSOMES}/{BIN_SIZE}/combined.merged.vcf.reformatted.single.event.only.melt.dt.txt"
    shell:
        """
        echo CHROM,POS,ID,REF,ALT,SUBSET,sample,AC,AN,AF > {output.combined_merged_vcf_reformatted_single_event_only_melt_dt_filename}
        for file in {input.merged_vcf_reformatted_single_event_only_melt_dt_filename_collection:q}
        do
            tail -n +2 $file >> {output.combined_merged_vcf_reformatted_single_event_only_melt_dt_filename}
        done
        """

def __S18_1_step02_write_annotation_dt_filename_collection_to_file_and_return_filename(wildcards):
    temp_result_directory = "result/S18_1__combine_annotations/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/snpEff/" + wildcards.CHROMOSOMES + "/" + wildcards.BIN_SIZE + "/"
    temp_annotation_dt_filename_list_filename = temp_result_directory + "/temp_annotation_dt_filename_list.txt"
    os.makedirs(temp_result_directory, exist_ok=True)
    temp_annotation_dt_filenames_collection = ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/snpEff/merged.vcf.reformatted.single.event.only.melt.with.snpEff.ANN.split.annotation.dt.txt" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    with open(temp_annotation_dt_filename_list_filename, "w") as temp_filehandle:
        print("writing...")
        temp_filehandle.writelines([temp_annotation_dt_filename + "\n" for temp_annotation_dt_filename in temp_annotation_dt_filenames_collection])
    return temp_annotation_dt_filename_list_filename


rule S18_1__combine_annotations____step02__combine_merged_vcf_reformatted_single_event_only_melt_with_snpEff_ANN_split_annotation_dt_filename:
    input:
        S17_2_step03_flags_collection=lambda wildcards: ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/snpEff/finished.step03__compute_tables_and_stats_based_on_snpEff_annotation" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    output:
        flag_S18_1_step02=touch("result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/{CHROMOSOMES}/{BIN_SIZE}/finished.step02__combine_merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_filename"),
        combined_editing_site_and_transcript_and_gene_dt_filename="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/snpEff/{COMPLEX_FILTER_PARAMETERS}/{CHROMOSOMES}/{BIN_SIZE}/combined.editing.site.and.transcript.and.gene.dt.txt"
    params:
        merged_vcf_reformatted_single_event_only_melt_with_snpEff_ANN_split_annotation_dt_filename_list_filename=lambda wildcards: __S18_1_step02_write_annotation_dt_filename_collection_to_file_and_return_filename(wildcards)
    shell:
        """
        echo CHROM,POS,ID,REF,ALT,SUBSET,sample,AC,AN,AF,Gene_Name,Gene_ID,Feature_ID,Transcript_BioType,event > {output.combined_editing_site_and_transcript_and_gene_dt_filename}
        cat {params.merged_vcf_reformatted_single_event_only_melt_with_snpEff_ANN_split_annotation_dt_filename_list_filename} | while read temp_annotation_filename
        do
            echo `date` Processing $temp_annotation_filename
            tail -n +2 $temp_annotation_filename | cut -d , -f 1,2,3,4,5,6,7,8,9,10,14,15,17,18,27 >> {output.combined_editing_site_and_transcript_and_gene_dt_filename}
        done
        """






def __S18_1_step02_patch01_write_merged_variant_only_snpEff_ANN_single_match_dt_filename_collection_to_file_and_return_filename(wildcards):
    temp_result_directory = "result/S18_1__combine_annotations/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/snpEff/" + wildcards.CHROMOSOMES + "/" + wildcards.BIN_SIZE + "/"
    temp_annotation_dt_filename_list_filename = temp_result_directory + "/temp_merged.variant.only.snpEff.ANN.single.match.dt.txt_list.txt"
    os.makedirs(temp_result_directory, exist_ok=True)
    temp_annotation_dt_filenames_collection = ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/snpEff/merged.variant.only.snpEff.ANN.single.match.dt.txt" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    with open(temp_annotation_dt_filename_list_filename, "w") as temp_filehandle:
        print("writing...")
        temp_filehandle.writelines([temp_annotation_dt_filename + "\n" for temp_annotation_dt_filename in temp_annotation_dt_filenames_collection])
    return temp_annotation_dt_filename_list_filename


rule S18_1__combine_annotations____step02__combine_merged_vcf_reformatted_single_event_only_melt_with_snpEff_ANN_split_annotation_dt_filename____patch01__get_full_annotation:
    input:
        S17_2_step03_flags_collection=lambda wildcards: ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/snpEff/finished.step03__compute_tables_and_stats_based_on_snpEff_annotation" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    output:
        flag_S18_1_step02=touch("result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/{CHROMOSOMES}/{BIN_SIZE}/finished.step02__combine_merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_filename____patch01__get_full_annotation"),
        combined_merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_gz_filename="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/{CHROMOSOMES}/{BIN_SIZE}/combined.merged.variant.only.snpEff.ANN.single.match.dt.txt.gz"
    params:
        merged_variant_only_snpEff_ANN_single_match_dt_filename_list_filename=lambda wildcards: __S18_1_step02_patch01_write_merged_variant_only_snpEff_ANN_single_match_dt_filename_collection_to_file_and_return_filename(wildcards),
        temp_directory="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/{CHROMOSOMES}/{BIN_SIZE}/temp_directory/"
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        rm -fr {params.temp_directory}/step02.patch.01.temp.result.txt
        mkdir -p {params.temp_directory}
        echo CHROM,POS,ID,REF,ALT,Allele,Annotation,Annotation_Impact,Gene_Name,Gene_ID,Feature_Type,Feature_ID,Transcript_BioType,Rank,HGVS.c,HGVS.p,cDNA.pos / cDNA.length,CDS.pos / CDS.length,AA.pos / AA.length,Distance,ERRORS / WARNINGS / INFO,event > {params.temp_directory}/step02.patch.01.temp.result.txt
        cat {params.merged_variant_only_snpEff_ANN_single_match_dt_filename_list_filename} | while read temp_annotation_filename
        do
            echo `date` Processing $temp_annotation_filename
            tail -n +2 $temp_annotation_filename >> {params.temp_directory}/step02.patch.01.temp.result.txt
        done
        pigz --processes {threads} {params.temp_directory}/step02.patch.01.temp.result.txt --stdout > {output.combined_merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_gz_filename}
        """








def __S18_1_step04_write_vcf_gz_filename_collection_to_file_and_return_filename(wildcards):
    temp_result_directory = "result/S18_1__combine_annotations/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/bcftools.isec.with.outer.vcf/" + wildcards.CHROMOSOMES + "/" + wildcards.BIN_SIZE + "/" + wildcards.OUTER_VCF_NAME + "/"
    temp_vcf_gz_filename_list_filename = temp_result_directory + "/temp_vcf_gz_filename_list.txt"
    os.makedirs(temp_result_directory, exist_ok=True)
    temp_vcf_gz_filenames_collection = ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/bcftools.isec.with.outer.vcf/" + wildcards.OUTER_VCF_NAME + "/" + temp_chromosome + "/merged.variant.only.bcftools.isec.outer.vcf.result.vcf.gz" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    with open(temp_vcf_gz_filename_list_filename, "w") as temp_filehandle:
        print("writing...")
        temp_filehandle.writelines([temp_vcf_gz_filename + "\n" for temp_vcf_gz_filename in temp_vcf_gz_filenames_collection])
    return temp_vcf_gz_filename_list_filename


rule S18_1__combine_annotations____step03__combine_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename:
    input:
        S17_2_step04_flags_collection=lambda wildcards: ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/bcftools.isec.with.outer.vcf/" +  wildcards.OUTER_VCF_NAME + "/" + temp_chromosome + "/" + "finished.step04__find_overlaps_with_known_DNA_variants" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    output:
        flag_S18_1_step02=touch("result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/bcftools.isec.with.outer.vcf/{CHROMOSOMES}/{BIN_SIZE}/{OUTER_VCF_NAME}/finished.step04__find_overlaps_with_known_DNA_variants"),
        combined_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/bcftools.isec.with.outer.vcf/{CHROMOSOMES}/{BIN_SIZE}/{OUTER_VCF_NAME}/combined.merged.variant.only.bcftools.isec.outer.vcf.result.vcf.gz"
    params:
        merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename_list_filename=lambda wildcards: __S18_1_step04_write_vcf_gz_filename_collection_to_file_and_return_filename(wildcards)
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        bcftools concat --threads {threads} --allow-overlaps -o {output.combined_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename} -Oz --file-list {params.merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename_list_filename}
        tabix {output.combined_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename}
        """



def __S18_1_step04_write_vcf_gz_filename_collection_to_file_and_return_filename____patch01__collapse_all(wildcards):
    temp_result_directory = "result/S18_1__combine_annotations/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/bcftools.isec.with.outer.vcf/" + wildcards.CHROMOSOMES + "/" + wildcards.BIN_SIZE + "/" + wildcards.OUTER_VCF_NAME + "/collapse_all/"
    temp_vcf_gz_filename_list_filename = temp_result_directory + "/temp_vcf_gz_filename_list.txt"
    os.makedirs(temp_result_directory, exist_ok=True)
    temp_vcf_gz_filenames_collection = ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/bcftools.isec.with.outer.vcf/" + wildcards.OUTER_VCF_NAME + "/" + temp_chromosome + "/collapse_all/merged.variant.only.bcftools.isec.outer.vcf.result.vcf.gz" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    with open(temp_vcf_gz_filename_list_filename, "w") as temp_filehandle:
        print("writing...")
        temp_filehandle.writelines([temp_vcf_gz_filename + "\n" for temp_vcf_gz_filename in temp_vcf_gz_filenames_collection])
    return temp_vcf_gz_filename_list_filename


rule S18_1__combine_annotations____step03__combine_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename____patch01__collapse_all:
    input:
        S17_2_step05_flags_collection=lambda wildcards: ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/bcftools.isec.with.outer.vcf/" +  wildcards.OUTER_VCF_NAME + "/" + temp_chromosome + "/collapse_all/" + "finished.step05__find_overlaps_with_known_DNA_variants_with_collapse_all" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    output:
        flag_S18_1_step02__patch01=touch("result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/bcftools.isec.with.outer.vcf/{CHROMOSOMES}/{BIN_SIZE}/{OUTER_VCF_NAME}/finished.step05__find_overlaps_with_known_DNA_variants_with_collapse_all"),
        combined_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/bcftools.isec.with.outer.vcf/{CHROMOSOMES}/{BIN_SIZE}/{OUTER_VCF_NAME}/collapse_all/combined.merged.variant.only.bcftools.isec.outer.vcf.result.vcf.gz"
    params:
        merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename_list_filename=lambda wildcards: __S18_1_step04_write_vcf_gz_filename_collection_to_file_and_return_filename____patch01__collapse_all(wildcards)
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        bcftools concat --threads {threads} --allow-overlaps -o {output.combined_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename} -Oz --file-list {params.merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename_list_filename}
        tabix {output.combined_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename}
        """

def __S19_1_step06_write_vcf_gz_filename_collection_to_file_and_return_filename_with_collapse_all_and_keep_self_vcf(wildcards):
    temp_result_directory = "result/S19_1__combine_annotations/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/bcftools.isec.with.outer.vcf/" + wildcards.CHROMOSOMES + "/" + wildcards.BIN_SIZE + "/" + wildcards.OUTER_VCF_NAME + "/collapse_all_and_keep_self_vcf/"
    temp_self_vcf_gz_filename_list_filename = temp_result_directory + "/temp_self_vcf_gz_filename_list.txt"
    os.makedirs(temp_result_directory, exist_ok=True)
    temp_self_vcf_gz_filenames_collection = ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/bcftools.isec.with.outer.vcf/" + wildcards.OUTER_VCF_NAME + "/" + temp_chromosome + "/collapse_all_and_keep_self_vcf/merged.variant.only.bcftools.isec.outer.vcf.result.self.vcf.gz" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    with open(temp_self_vcf_gz_filename_list_filename, "w") as temp_self_filehandle:
        print("writing...")
        temp_self_filehandle.writelines([temp_self_vcf_gz_filename + "\n" for temp_self_vcf_gz_filename in temp_self_vcf_gz_filenames_collection])
    return temp_self_vcf_gz_filename_list_filename


rule S19_1__combine_annotations____step06__combine_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename_with_collapse_all_and_keep_self_vcf:
    input:
        S17_2_step05_flags_collection=lambda wildcards: ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/bcftools.isec.with.outer.vcf/" +  wildcards.OUTER_VCF_NAME + "/" + temp_chromosome + "/collapse_all_and_keep_self_vcf/" + "finished.step06__find_overlaps_with_known_DNA_variants_with_collapse_all_and_keep_self_vcf" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    output:
        flag_S19_1_step06=touch("result/S19_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/bcftools.isec.with.outer.vcf/{CHROMOSOMES}/{BIN_SIZE}/{OUTER_VCF_NAME}/collapse_all_and_keep_self_vcf/finished.step06__combine_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_vcf_gz_filename_with_collapse_all_and_keep_self_vcf"),
        combined_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_self_vcf_gz_filename="result/S19_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/bcftools.isec.with.outer.vcf/{CHROMOSOMES}/{BIN_SIZE}/{OUTER_VCF_NAME}/collapse_all_and_keep_self_vcf/combined.merged.variant.only.bcftools.isec.outer.vcf.result.self.vcf.gz"
    params:
        merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_self_vcf_gz_filename_list_filename=lambda wildcards: __S19_1_step06_write_vcf_gz_filename_collection_to_file_and_return_filename_with_collapse_all_and_keep_self_vcf(wildcards)
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        bcftools concat --threads {threads} --allow-overlaps -o {output.combined_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_self_vcf_gz_filename} -Oz --file-list {params.merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_self_vcf_gz_filename_list_filename}
        tabix {output.combined_merged_variant_only_vcf_gz_bcftools_isec_outer_vcf_result_self_vcf_gz_filename}
        """



def __S18_1_step05_write_merged_variant_only_snpEff_event_summary_dt_filename_collection_to_file_and_return_filename(wildcards):
    temp_result_directory = "result/S18_1__combine_annotations/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/snpEff/" + wildcards.CHROMOSOMES + "/" + wildcards.BIN_SIZE + "/"
    temp_annotation_dt_filename_list_filename = temp_result_directory + "/temp_merged_variant_only_snpEff_event_summary_dt_filename_list.txt"
    os.makedirs(temp_result_directory, exist_ok=True)
    temp_annotation_dt_filenames_collection = ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/snpEff/merged.variant.only.snpEff.event.summary.dt.txt" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    with open(temp_annotation_dt_filename_list_filename, "w") as temp_filehandle:
        print("writing...")
        temp_filehandle.writelines([temp_annotation_dt_filename + "\n" for temp_annotation_dt_filename in temp_annotation_dt_filenames_collection])
    return temp_annotation_dt_filename_list_filename


rule S18_1__combine_annotations____step05__combine_merged_variant_only_snpEff_event_summary_dt_filename:
    input:
        S17_2_step03_flags_collection=lambda wildcards: ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/snpEff/finished.step03__compute_tables_and_stats_based_on_snpEff_annotation" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    output:
        flag_S18_1_step02=touch("result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/{CHROMOSOMES}/{BIN_SIZE}/finished.step05__combine_merged_variant_only_snpEff_event_summary_dt_filename"),
        combined_merged_variant_only_snpEff_event_summary_dt_gz_filename="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/{CHROMOSOMES}/{BIN_SIZE}/combined.merged.variant.only.snpEff.event.summary.dt.txt.gz"
    params:
        merged_variant_only_snpEff_event_summary_dt_filename_list_filename=lambda wildcards: __S18_1_step05_write_merged_variant_only_snpEff_event_summary_dt_filename_collection_to_file_and_return_filename(wildcards),
        temp_directory="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/{CHROMOSOMES}/{BIN_SIZE}/temp_directory/"
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        rm -fr {params.temp_directory}/step05.temp.result.txt
        mkdir -p {params.temp_directory}
        echo CHROM,POS,ID,REF,ALT,event.summary > {params.temp_directory}/step05.temp.result.txt
        cat {params.merged_variant_only_snpEff_event_summary_dt_filename_list_filename} | while read temp_annotation_filename
        do
            echo `date` Processing $temp_annotation_filename
            tail -n +2 $temp_annotation_filename >> {params.temp_directory}/step05.temp.result.txt
        done
        pigz --processes {threads} --stdout {params.temp_directory}/step05.temp.result.txt > {output.combined_merged_variant_only_snpEff_event_summary_dt_gz_filename}
        """




def __S18_1_step06_write_merged_vcf_reformatted_single_event_only_melt_dt_filename_collection_to_file_and_return_filename(wildcards):
    temp_result_directory = "result/S18_1__combine_annotations/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/snpEff/" + wildcards.CHROMOSOMES + "/" + wildcards.BIN_SIZE + "/"
    temp_annotation_dt_filename_list_filename = temp_result_directory + "/temp_merged_vcf_reformatted_single_event_only_melt_dt_filename_list.txt"
    os.makedirs(temp_result_directory, exist_ok=True)
    temp_annotation_dt_filenames_collection = ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/snpEff/merged.vcf.reformatted.single.event.only.melt.dt.txt" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    with open(temp_annotation_dt_filename_list_filename, "w") as temp_filehandle:
        print("writing...")
        temp_filehandle.writelines([temp_annotation_dt_filename + "\n" for temp_annotation_dt_filename in temp_annotation_dt_filenames_collection])
    return temp_annotation_dt_filename_list_filename


rule S18_1__combine_annotations____step06__combine_merged_vcf_reformatted_single_event_only_melt_dt_filename:
    input:
        S17_2_step03_flags_collection=lambda wildcards: ["result/S17_2__annotate_merged_sites_from_a_dataset_collection_for_a_single_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_chromosome + ":" + str(temp_bin_start_closed) + "-" + str(temp_bin_end_closed) + "/snpEff/finished.step03__compute_tables_and_stats_based_on_snpEff_annotation" for temp_chromosome, temp_bin_start_closed, temp_bin_end_closed in __get_uniform_regions(wildcards.CHROMOSOMES, wildcards.BIN_SIZE)]
    output:
        flag_S18_1_step02=touch("result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/snpEff/{CHROMOSOMES}/{BIN_SIZE}/finished.step06__combine_merged_vcf_reformatted_single_event_only_melt_dt_filename"),
        combined_merged_vcf_reformatted_single_event_only_melt_dt_filename="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/snpEff/{COMPLEX_FILTER_PARAMETERS}/{CHROMOSOMES}/{BIN_SIZE}/combined.merged.vcf.reformatted.single.event.only.melt.dt.txt"
    params:
        merged_vcf_reformatted_single_event_only_melt_dt_filename_list_filename=lambda wildcards: __S18_1_step06_write_merged_vcf_reformatted_single_event_only_melt_dt_filename_collection_to_file_and_return_filename(wildcards)
    shell:
        """
        echo CHROM,POS,ID,REF,ALT,event.summary > {output.combined_merged_variant_only_snpEff_event_summary_dt_filename}
        cat {params.merged_vcf_reformatted_single_event_only_melt_dt_filename_list_filename} | while read temp_annotation_filename
        do
            echo `date` Processing $temp_annotation_filename
            tail -n +2 $temp_annotation_filename >> {output.combined_merged_vcf_reformatted_single_event_only_melt_dt_filename}
        done
        """



rule S06_1__get_expression_level__step01__align_sample:
    input:
        flag_S01_2__step01="result/S01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/finished.step01__trim_RNA_Seq_by_fastp",
        flag_s04_1="result/s05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/finished"
    output:
        flag_S06_1__step01=touch("result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01__align_sample")
    params:
        contigs_index_prefix="result/s05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/prefix",
        trim_result_directory="result/S01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/",
        alignment_result_directory="result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/"
    threads:
        config['threads_aligning']
    shell:
        """
        type_end=`echo {wildcards.TYPE} | cut -f 1 -d '-'`
        rm -fr "{params.alignment_result_directory}/"*
        mkdir -p "{params.alignment_result_directory}/"
        if [[ "{wildcards.ALIGNER}" == 'STAR-expression' ]]; then
            aligner_parameters="{wildcards.ALIGNER_PARAMETERS}"
            star_parameters="-h"
            if [[ $aligner_parameters == 'default' ]]; then
                star_parameters=""
            fi

            if [[ $type_end == single ]]; then
                STAR --runThreadN {threads} --genomeDir {params.contigs_index_prefix} --readFilesIn "{params.trim_result_directory}/r.fastp.fastq.gz"   --readFilesCommand zcat --outFileNamePrefix {params.alignment_result_directory}/out_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --chimSegmentMin 20 --outSAMstrandField intronMotif --chimOutType SeparateSAMold $star_parameters
                ln -s out_Aligned.sortedByCoord.out.bam {params.alignment_result_directory}/alignment.bam
            elif [[ $type_end == paired ]]; then
                STAR --runThreadN {threads} --genomeDir {params.contigs_index_prefix} --readFilesIn "{params.trim_result_directory}/r1.fastp.fastq.gz" "{params.trim_result_directory}/r2.fastp.fastq.gz" --readFilesCommand zcat --outFileNamePrefix {params.alignment_result_directory}/out_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --chimSegmentMin 20 --outSAMstrandField intronMotif --chimOutType SeparateSAMold $star_parameters
                ln -s out_Aligned.sortedByCoord.out.bam {params.alignment_result_directory}/alignment.bam
            else
                echo 'Unsupported read type {wildcards.TYPE}' && exit 1
            fi
        else
            echo 'Unsupported ALIGNER {wildcards.ALIGNER}' && exit 1
        fi
        """


rule S06_1__get_expression_level__step01_2__merge_bams_of_runs_for_the_same_sample:
    input:
        flags_collection_S06_1__step01 = lambda wildcards: ["result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/" + temp_run + "/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01__align_sample" for temp_run in __get_RUNs_by_DATASET_and_SAMPLE(wildcards.DATASET, wildcards.SAMPLE)]
    output:
        flag_S06_1__step01_2=touch("result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01_2__merge_bams_of_runs_for_the_same_sample"),
        alignment_merged_filename="result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/alignment.merged.bam"
    params:
        alignment_bams_prefix = "result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/",
        alignment_bams_suffix = "/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/alignment.bam",
        alignment_bams_runs = lambda wildcards: __get_RUNs_by_DATASET_and_SAMPLE(wildcards.DATASET, wildcards.SAMPLE),
        number_of_alignments_to_merge = lambda wildcards: len(__get_RUNs_by_DATASET_and_SAMPLE(wildcards.DATASET, wildcards.SAMPLE)),
        merge_result_directory="result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/",
        picard_Xmx=config["picard_Xmx"]
    threads:
        config['threads_merging_bams']
    shell:
        """
        if [[ {params.number_of_alignments_to_merge} > 1 ]]; then
            final_INPUT_literal=''
            for item in {params.alignment_bams_runs:q}
            do
                final_INPUT_literal=$final_INPUT_literal"INPUT="{params.alignment_bams_prefix}$item{params.alignment_bams_suffix}" "
            done
            picard MergeSamFiles {params.picard_Xmx} $final_INPUT_literal OUTPUT={output.alignment_merged_filename} SO=coordinate VALIDATION_STRINGENCY=SILENT
        else ## == 1
            ln -s `readlink -f {params.alignment_bams_prefix}{params.alignment_bams_runs}{params.alignment_bams_suffix}` {output.alignment_merged_filename}
        fi
        """




rule BS06_1__get_expression_level__step02__call_expression:
    input:
        lambda wildcards: ["result/S06_1__get_expression_level/" + temp_row.TYPE + '/' + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/" + str(temp_row.INDEXER_PARAMETERS) + "/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/finished.step02__call_expression" for temp_row in __get_parameters_collection_for_expression_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_BS06_1__step02=touch("result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/finished.step02__call_expression")

rule BS06_1__get_expression_level__step03__get_expression_matrix_by_ballgown:
    input:
        flag_BS06_1__step02="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/finished.step02__call_expression"
    output:
        flag_S06_1__step03=touch("result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/finished.step03__get_expression_matrix_by_ballgown"),
        combined_texpr_FPKM_matrix_filename="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/combined.texpr.FPKM.matrix.txt",
        combined_gexpr_FPKM_matrix_filename="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/combined.gexpr.FPKM.matrix.txt"
    params:
        temp_linking_directory="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/temp_linking/",
        ballgown_directory_filenames_collection = lambda wildcards: list(dict.fromkeys(["result/S06_1__get_expression_level/" + temp_row.TYPE + '/' + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/__merged__/" + wildcards.FASTP_PARAMETERS + '/' + wildcards.CONTIGS_FASTA_FILENAME + '/' + wildcards.GENCODE_VERSION + '/' + wildcards.INDEXER + '/' + str(temp_row.INDEXER_PARAMETERS) + '/' + wildcards.ALIGNER + '/' + wildcards.ALIGNER_PARAMETERS + '/' + wildcards.EXPRESSION_CALLER + '/' + wildcards.EXPRESSION_CALLER_PARAMETERS + '/' for temp_row in __get_parameters_collection_for_expression_calling(wildcards.DATASET_COLLECTION_NAME)])) ## must remove GSM-level duplicates (here with order preserved); otherwise the ballgown will throw an error
    threads:
        config["threads_auxiliary_processing"]
    script:
        "./scripts/BS06_1__get_expression_level/step03__get_expression_matrix/run.R"





rule BS06_1__get_expression_level__step04__get_count_matrix:
    input:
        flag_BS06_1__step02="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/finished.step02__call_expression"
    output:
        flag_S06_1__step04=touch("result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/finished.step04__get_count_matrix_by_ballgown"),
        combined_texpr_cov_matrix_filename="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/combined.texpr.count.matrix.txt",
        combined_gexpr_cov_matrix_filename="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/combined.gexpr.count.matrix.txt"
    params:
        temp_linking_directory="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/step04__temp_linking/",
        ballgown_directory_filenames_collection = lambda wildcards: list(dict.fromkeys(["result/S06_1__get_expression_level/" + temp_row.TYPE + '/' + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/__merged__/" + wildcards.FASTP_PARAMETERS + '/' + wildcards.CONTIGS_FASTA_FILENAME + '/' + wildcards.GENCODE_VERSION + '/' + wildcards.INDEXER + '/' + str(temp_row.INDEXER_PARAMETERS) + '/' + wildcards.ALIGNER + '/' + wildcards.ALIGNER_PARAMETERS + '/' + wildcards.EXPRESSION_CALLER + '/' + wildcards.EXPRESSION_CALLER_PARAMETERS + '/' for temp_row in __get_parameters_collection_for_expression_calling(wildcards.DATASET_COLLECTION_NAME)])) ## must remove GSM-level duplicates (here with order preserved); otherwise the ballgown will throw an error
    threads:
        config["threads_auxiliary_processing"]
    script:
        "./scripts/BS06_1__get_expression_level/step04__get_count_matrix/run.R"






rule S20_1__extract_GSE_table_by_GEOmetadb:
    input:
        GEOmetadb_sqlite_filename="external/NCBI.GEOmetadb/GEOmetadb.sqlite"
    output:
        flag_S20_1=touch("result/S20_1__extract_GSE_table_by_GEOmetadb/{GSE}/finished"),
        gse_and_gsm_dt_filename="result/S20_1__extract_GSE_table_by_GEOmetadb/{GSE}/gse.and.gsm.dt.txt"
    script:
        "./scripts/S20_1__extract_GSE_table_by_GEOmetadb/run.R"


def __get_all_gse_from_dataset_phenotype_collection(dataset_phenotype_collection_name, prefix_for_dataset_phenotype_collection_name_directory="external/DATASET_PHENOTYPE_COLLECTION_NAME_DIRECTORY/"):
    temp_dataset_DataFrame = pandas.read_csv( prefix_for_dataset_phenotype_collection_name_directory + "/" +  dataset_phenotype_collection_name, sep=",")
    print(temp_dataset_DataFrame)
    temp_GSE_collection = [ temp_row["STUDY"] for temp_rowindex, temp_row in temp_dataset_DataFrame.iterrows()]
    return temp_GSE_collection


rule S21_1__merge_phenotype_tables:
    input:
        phenotype_collection_filename="external/DATASET_PHENOTYPE_COLLECTION_NAME_DIRECTORY/{DATASET_PHENOTYPE_COLLECTION_NAME}",
        gse_and_gsm_dt_filename_collection= lambda wildcards: ["result/S20_1__extract_GSE_table_by_GEOmetadb/" + temp_GSE + "/gse.and.gsm.dt.txt" for temp_GSE in __get_all_gse_from_dataset_phenotype_collection(wildcards.DATASET_PHENOTYPE_COLLECTION_NAME)],
        manuscript_table_for_processing_xlsx_filename="manuscript/table_for_processing.xlsx"
    output:
        flag_S21_1=touch("result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished"),
        phenotype_output_dt_filename="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/phenotype.output.dt.txt",
        phenotype_output_at_gsm_level_dt_filename="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/phenotype.output.at.gsm.level.dt.txt"
    threads:
        config["threads_auxiliary_processing"]
    script:
        "./scripts/S21_1__merge_phenotype_tables/run.R"






rule S41_1__check_variant_converage_of_merged_bam____step01__generate_bed_from_vcf:
    input:
        flag_S16_1="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished",
        merged_variant_only_vcf_gz_filename="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/merged.variant.only.vcf.gz"
    output:
        flag_S41_1_step01=touch("result/S41_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step01__generate_bed_from_vcf"),
        merged_variant_only_bed_filename="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/merged.variant.only.bed"
    threads:
        1
    shell:
        """
        bcftools view {input.merged_variant_only_vcf_gz_filename} | vcf2bed > {output.merged_variant_only_bed_filename}
        """



rule S41_1__check_variant_converage_of_merged_bam____step02__check_a_single_bam_coverage_of_the_bed_generated_from_vcf:
    input:
        flag_S41_1_step01="result/S41_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step01__generate_bed_from_vcf",
        merged_variant_only_bed_filename="result/S41_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/merged.variant.only.bed",
        flag_S15_1__step01_2="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01_2__merge_bams_of_runs_for_the_same_sample",
        alignment_merged_filename="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/alignment.merged.bam"
    output:
        flag_S41_1_step02=touch("result/S41_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished.step02__check_a_single_bam_coverage_of_the_bed_generated_from_vcf"),
        alignment_merged_depth_on_merged_variant_only_bed_txt_gz_filename="result/S41_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.merged.variant.only.bed.txt.gz"
    threads:
        1
    shell:
        """
        samtools depth -a -d 0 -b {input.merged_variant_only_bed_filename} {input.alignment_merged_filename} | gzip > {output.alignment_merged_depth_on_merged_variant_only_bed_txt_gz_filename}
        """



rule B41_1__check_variant_converage_of_merged_bam____step02__check_a_single_bam_coverage_of_the_bed_generated_from_vcf:
    input:
        lambda wildcards: ["result/S41_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/finished.step02__check_a_single_bam_coverage_of_the_bed_generated_from_vcf" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME_FOR_READ_COVERAGE_CHECK)]
    output:
        flag_B41_1__step02=touch("result/B41_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/{DATASET_COLLECTION_NAME_FOR_READ_COVERAGE_CHECK}/finished.step02__check_a_single_bam_coverage_of_the_bed_generated_from_vcf")
