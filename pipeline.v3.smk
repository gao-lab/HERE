import pandas
import re

for temp_key in ['threads_qc', 'threads_indexing', 'threads_trimming', 'threads_aligning', 'threads_calling_variants', 'threads_merging_bams', 'threads_calling_expression', 'threads_merging_vcfs', 'threads_subset_bam']:
    if temp_key not in config.keys():
        config[temp_key] = 20

for temp_key in ['threads_calling_expression']:
    if temp_key not in config.keys():
        config[temp_key] = 4

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

rule S01_1__check_RNA_Seq____step01__check_RNA_Seq_by_fastqc:
    input:
        input_read_filenames=__get_raw_fastq_path_by_TYPE_and_DATASET_and_SAMPLE_and_RUN
    output:
        flag_S01_2__step01=touch("result/S01_1__check_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTQC_PARAMETERS}/finished.step01__check_RNA_Seq_by_fastqc")
    params:
        raw_fastq_parent_path="external/RNA-Seq-with-Run/{DATASET}/{SAMPLE}/{RUN}/RNA/",
        fastqc_result_directory="result/S01_1__check_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTQC_PARAMETERS}/fastqc_result/"
    threads:
        config['threads_qc']
    shell:
        """
        type_end=`echo {wildcards.TYPE} | cut -f 1 -d '-'`
        rm -fr "{params.fastqc_result_directory}"/* && mkdir -p "{params.fastqc_result_directory}"
        fastqc_parameters_literal="{wildcards.FASTQC_PARAMETERS}"
        fastqc_parameters=""
        if [[ $fastqc_parameters_literal == 'default' ]]; then
            fastqc_parameters=""
        else
            echo "Unsupported fastqc parameters $fastqc_parameters" && exit 1
        fi

        if [[ $type_end == 'paired' ]]; then
            mkdir -p "{params.fastqc_result_directory}"/r1/
            fastqc $fastqc_parameters --outdir "{params.fastqc_result_directory}"/r1/ --threads {threads} {params.raw_fastq_parent_path}/r1.fastq.gz
            mkdir -p "{params.fastqc_result_directory}"/r2/
            fastqc $fastqc_parameters --outdir "{params.fastqc_result_directory}"/r2/ --threads {threads} {params.raw_fastq_parent_path}/r2.fastq.gz
        elif [[ $type_end == 'single' ]]; then
            mkdir -p "{params.fastqc_result_directory}"/r/
            fastqc $fastqc_parameters --outdir "{params.fastqc_result_directory}"/r/ --threads {threads} {params.raw_fastq_parent_path}/r.fastq.gz
        else
            echo "Unsupported type: {wildcards.TYPE}" && exit 1
        fi

        """


rule B01_1__check_RNA_Seq____step01__check_RNA_Seq_by_fastqc:
    input:
        lambda wildcards: ["result/S01_1__check_RNA_Seq/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/" + temp_RUN + "/{FASTQC_PARAMETERS}/finished.step01__check_RNA_Seq_by_fastqc" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME) for temp_RUN in __get_RUNs_by_DATASET_and_SAMPLE(temp_row.DATASET_NAME, temp_row.SAMPLE_NAME)]
    output:
        flag_B01_1__step01=touch("result/B01_1__check_RNA_Seq/{DATASET_COLLECTION_NAME}/__merged__/{FASTQC_PARAMETERS}/finished.step01__check_RNA_Seq_by_fastqc")

rule S01_1__check_RNA_Seq____step02__extract_overlapping_sequences_and_adapters:
    input:
        flag_B01_1__step01="result/B01_1__check_RNA_Seq/{DATASET_COLLECTION_NAME}/__merged__/{FASTQC_PARAMETERS}/finished.step01__check_RNA_Seq_by_fastqc"
    output:
        flag_S01_1__step02=touch("result/S01_1__check_RNA_Seq/{DATASET_COLLECTION_NAME}/__merged__/{FASTQC_PARAMETERS}/finished.step02__extract_overlapping_sequences_and_adapters"),
        adapter_content_dt_txt_filename="result/S01_1__check_RNA_Seq/{DATASET_COLLECTION_NAME}/__merged__/{FASTQC_PARAMETERS}/adapter.content.dt.txt",
        overrepresented_sequences_dt_txt_filename="result/S01_1__check_RNA_Seq/{DATASET_COLLECTION_NAME}/__merged__/{FASTQC_PARAMETERS}/overrepresented.sequences.dt.txt"
    params:
        fastqc_result_directory_collection=lambda wildcards: ["result/S01_1__check_RNA_Seq/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/" + temp_RUN + "/" + wildcards.FASTQC_PARAMETERS + "/fastqc_result/" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME) for temp_RUN in __get_RUNs_by_DATASET_and_SAMPLE(temp_row.DATASET_NAME, temp_row.SAMPLE_NAME)]
    script:
        "./scripts/S01_1__check_RNA_Seq/step02__extract_overlapping_sequences_and_adapters/run.R"


rule S01_2__trim_RNA_Seq____step01__trim_RNA_Seq:
    input:
        input_read_filenames=__get_raw_fastq_path_by_TYPE_and_DATASET_and_SAMPLE_and_RUN
    output:
        flag_S01_2__step01=touch("result/S01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/finished.step01__trim_RNA_Seq")
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
            elif [[ "{wildcards.FASTP_PARAMETERS}" == 'base-quality-no-smaller-than-25' ]]; then ##the quality filter did not effect due to --disable_quality_filtering
                fastp --disable_adapter_trimming --disable_trim_poly_g --disable_quality_filtering --average_qual 25 --disable_length_filtering  -i {params.raw_fastq_parent_path}/r1.fastq.gz -I {params.raw_fastq_parent_path}/r2.fastq.gz -o "{params.trim_result_directory}/r1.fastp.fastq.gz" -O "{params.trim_result_directory}/r2.fastp.fastq.gz" -w {threads} {wildcards.FASTP_PARAMETERS}
            elif [[ "{wildcards.FASTP_PARAMETERS}" == 'auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp' ]]; then
              trim_galore --quality 0 --paired --output_dir {params.trim_result_directory}/ --basename trim_galore --cores {threads} {params.raw_fastq_parent_path}/r1.fastq.gz {params.raw_fastq_parent_path}/r2.fastq.gz
              fastp --disable_adapter_trimming --average_qual 25  -i {params.trim_result_directory}/trim_galore_val_1.fq.gz -I {params.trim_result_directory}/trim_galore_val_2.fq.gz -o "{params.trim_result_directory}/r1.fastp.fastq.gz" -O "{params.trim_result_directory}/r2.fastp.fastq.gz" -w {threads}
            elif [[ "{wildcards.FASTP_PARAMETERS}" == 'DNTRSeq-DNA-trimming' ]]; then
                rm -fr {params.trim_result_directory}/* && mkdir -p {params.trim_result_directory}/
                trim_galore --paired --nextera -q 20 -e 0.1 --basename temp -o "{params.trim_result_directory}/" {params.raw_fastq_parent_path}/r1.fastq.gz {params.raw_fastq_parent_path}/r2.fastq.gz
                ln -s -r "{params.trim_result_directory}/"temp_val_1.fq.gz "{params.trim_result_directory}/"r1.fastp.fastq.gz
                ln -s -r "{params.trim_result_directory}/"temp_val_2.fq.gz "{params.trim_result_directory}/"r2.fastp.fastq.gz
            elif [[ "{wildcards.FASTP_PARAMETERS}" == 'DNTRSeq-RNA-trimming' ]]; then
                rm -fr {params.trim_result_directory}/* && mkdir -p {params.trim_result_directory}/temp
                zcat {params.raw_fastq_parent_path}/r1.fastq.gz > "{params.trim_result_directory}/temp/raw.r1.fastq"
                zcat {params.raw_fastq_parent_path}/r2.fastq.gz > "{params.trim_result_directory}/temp/raw.r2.fastq"
                prinseq-lite.pl -fastq "{params.trim_result_directory}/temp/raw.r1.fastq" -fastq2 "{params.trim_result_directory}/temp/raw.r2.fastq" -out_format 3 -trim_left 3 -min_len 15 -trim_qual_right 25 -lc_method entropy -lc_threshold 65 -out_good {params.trim_result_directory}/temp/prinseqRun1 -out_bad null
                prinseq-lite.pl -fastq {params.trim_result_directory}/temp/prinseqRun1_1.fastq -fastq2 {params.trim_result_directory}/temp/prinseqRun1_2.fastq -out_format 3 -min_len 15 -out_good {params.trim_result_directory}/temp/prinseqRun2 -out_bad null
                trim_galore --paired -a 'CTGTCTCTTATACACATCT' --length 15 --stringency 1 --fastqc_args "-j java --outdir={params.trim_result_directory}/fastqc_output --extract" --basename temp -o "{params.trim_result_directory}/" {params.trim_result_directory}/temp/prinseqRun2_1.fastq {params.trim_result_directory}/temp/prinseqRun2_2.fastq
                gzip -c "{params.trim_result_directory}/"temp_val_1.fq > "{params.trim_result_directory}/"r1.fastp.fastq.gz
                gzip -c "{params.trim_result_directory}/"temp_val_2.fq > "{params.trim_result_directory}/"r2.fastp.fastq.gz
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
            elif [[ "{wildcards.FASTP_PARAMETERS}" == 'auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp' ]]; then
              trim_galore --quality 0 --output_dir {params.trim_result_directory}/ --basename trim_galore --cores {threads} {params.raw_fastq_parent_path}/r.fastq.gz
              fastp --disable_adapter_trimming --average_qual 25  -i {params.trim_result_directory}/trim_galore_trimmed.fq.gz -o "{params.trim_result_directory}/r.fastp.fastq.gz"
            else
                fastp -i {params.raw_fastq_parent_path}/r.fastq.gz -o "{params.trim_result_directory}/r.fastp.fastq.gz" -w {threads} {wildcards.FASTP_PARAMETERS}
            fi
        else
            echo "Unsupported type: {wildcards.TYPE}" && exit 1
        fi

        """

rule s01_3__build_kallisto_index_for_check_strandedness:
    input:
        gencode_v32_transcripts_fa="external/contigs/gencode.v32.transcripts.fa"
    output:
        flag=touch("result/s01_3__build_kallisto_index_for_check_strandedness/finished")
    params:
        output_directory="result/s01_3__build_kallisto_index_for_check_strandedness"
    shell:
        """
        mkdir -p {params.output_directory}
        sed 's/[|]/ /g' {input.gencode_v32_transcripts_fa} > {params.output_directory}/gencode_v32_transcripts_renamed.fasta
        kallisto index -i {params.output_directory}/kallisto_index {params.output_directory}/gencode_v32_transcripts_renamed.fasta
        """

rule S01_3__check_strandedness:
    input:
        input_read_filenames=__get_raw_fastq_path_by_TYPE_and_DATASET_and_SAMPLE_and_RUN,
        gencode_v32_gtf="external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf",
        s01_3_flag="result/s01_3__build_kallisto_index_for_check_strandedness/finished"
    output:
        flag_S01_3=touch("result/S01_3__check_strandedness/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/finished")
    params:
        raw_fastq_parent_path="external/RNA-Seq-with-Run/{DATASET}/{SAMPLE}/{RUN}/RNA/",
        output_directory="result/S01_3__check_strandedness/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/",
        s01_3_renamed_fasta="result/s01_3__build_kallisto_index_for_check_strandedness/gencode_v32_transcripts_renamed.fasta",
        s01_3_kallisto_index="result/s01_3__build_kallisto_index_for_check_strandedness/kallisto_index"
    threads:
        1
    shell:
        """
        type_end=`echo {wildcards.TYPE} | cut -f 1 -d '-'`
        rm -fr "{params.output_directory}"/* && mkdir -p "{params.output_directory}"
        ln -s -r {input.gencode_v32_gtf} {params.output_directory}/gencode.v32.gtf
        ln -s -r {params.s01_3_renamed_fasta} {params.output_directory}/gencode_v32_transcripts_renamed.fasta
        ln -s -r {params.s01_3_kallisto_index} {params.output_directory}/kallisto_index
        ##
        if [[ $type_end == 'paired' ]]; then
            ln -s -r {params.raw_fastq_parent_path}/r1.fastq.gz {params.output_directory}/r1.fastq.gz
            ln -s -r {params.raw_fastq_parent_path}/r2.fastq.gz {params.output_directory}/r2.fastq.gz
            cd {params.output_directory}
            check_strandedness --gtf ./gencode.v32.gtf --transcripts ./gencode_v32_transcripts_renamed.fasta --kallisto_index ./kallisto_index --reads_1 ./r1.fastq.gz --reads_2 ./r2.fastq.gz
            ln -s ./stranded_test_r1/strandedness_check.txt .
        elif [[ $type_end == 'single' ]]; then
            ln -s -r {params.raw_fastq_parent_path}/r.fastq.gz {params.output_directory}/r.fastq.gz
            cd {params.output_directory}
            echo "single-end" > strandedness_check.txt 
        else
            echo "Unsupported type: {wildcards.TYPE}" && exit 1
        fi
        """
        
rule B01_3__check_strandedness:
    input:
        lambda wildcards: ["result/S01_3__check_strandedness/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/" + temp_RUN + "/finished" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME) for temp_RUN in __get_RUNs_by_DATASET_and_SAMPLE(temp_row.DATASET_NAME, temp_row.SAMPLE_NAME)]
    output:
        flag_B01_3__step01=touch("result/B01_3__check_strandedness/{DATASET_COLLECTION_NAME}/finished")


## 200831 updated version: replace 'bcftoolsxxx' with 'variant callers xxx'
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
        contigs_index_prefix="result/s05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/prefix",
        salmon_hg38_decoy_contig_fasta_filename="external/contigs/hg38.fa"
    threads:
        config['threads_indexing']
    shell:
        """
        rm -fr "{params.contigs_index_directory}"/*
        mkdir -p "{params.contigs_index_directory}"
        if [[ "{wildcards.INDEXER}" == 'bwa-index-default' ]]; then
            bwa index -p "{params.contigs_index_prefix}" "{input.contigs_fasta_filename}"
        elif [[ "{wildcards.INDEXER}" == 'bwa-index-10.1038_nmeth.2330' ]]; then
            cat "{params.reference_UCSC_Table_Browser_knownGene_filename}" | awk '{{print "bin\t"$0}}'  > "{params.contigs_index_directory}"/knownGene_with_dummy_bin_prepended
            perl ./scripts/s05_1__index_contig_with_annotation/getSplicingJunctions_modified.pl "{input.contigs_fasta_filename}" "{params.contigs_index_directory}"/knownGene_with_dummy_bin_prepended  {wildcards.INDEXER_PARAMETERS}  "{params.contigs_index_directory}"/splice_junctions.txt && cat "{input.contigs_fasta_filename}" "{params.contigs_index_directory}"/splice_junctions.txt > "{params.contigs_index_directory}"/combined.fasta
            bwa index -p "{params.contigs_index_prefix}" "{params.contigs_index_directory}"/combined.fasta
        elif [[ "{wildcards.INDEXER}" == 'STAR-expression' ]]; then
            rm -fr {params.contigs_index_prefix}/temp
            rm -fr {params.contigs_index_prefix}
            mkdir -p {params.contigs_index_prefix}
            STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.contigs_index_prefix} --genomeFastaFiles {input.contigs_fasta_filename} --outTmpDir {params.contigs_index_prefix}/temp --sjdbGTFfile {params.reference_GTF_filename} --sjdbOverhang {wildcards.INDEXER_PARAMETERS}
        elif [[ "{wildcards.INDEXER}" == 'salmon-gencode' ]]; then
            rm -fr {params.contigs_index_prefix}/*
            mkdir -p {params.contigs_index_prefix}
            cat {input.contigs_fasta_filename} {params.salmon_hg38_decoy_contig_fasta_filename} > {params.contigs_index_prefix}/transcripts.and.decoys.fasta
            grep ">" {params.salmon_hg38_decoy_contig_fasta_filename} | tr -d '>' > {params.contigs_index_prefix}/decoys.ids
            salmon index --transcripts {params.contigs_index_prefix}/transcripts.and.decoys.fasta --gencode --index {params.contigs_index_prefix} --threads {threads} --decoys {params.contigs_index_prefix}/decoys.ids
        else
            echo "Indexer {wildcards.INDEXER} not supported." && exit 1;
        fi
        """


rule S15_1__get_sample_RNA_editing_sites_v3__step01__align_sample:
    input:
        flag_S01_2__step01="result/S01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/finished.step01__trim_RNA_Seq",
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
        fastp_parameters="{wildcards.FASTP_PARAMETERS}"
        if [[ $fastp_parameters == "DNTRSeq-DNA-trimming" ]]; then
            samtools view -Sq 20 -h -b -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.bam" "{params.calling_result_directory}/alignment.sorted.withRG.dedup.bam"
            picard SortSam  {params.picard_Xmx}  INPUT="{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.bam" OUTPUT="{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=YES
            rm -f "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
            ln -s -r "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.bam" "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
            samtools index -@1 "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
        else
            samtools view -h "{params.calling_result_directory}/alignment.sorted.withRG.dedup.bam" > "{params.calling_result_directory}/alignment.sorted.withRG.dedup.sam"
            java {params.convertCoordinates_Xmx} -classpath tools/convertCoordinates_classpath convertCoordinates < "{params.calling_result_directory}/alignment.sorted.withRG.dedup.sam" | samtools view -b -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bam"
            rm -fr "{params.calling_result_directory}/alignment.sorted.withRG.dedup.sam"
            samtools view -Sq 20 -h "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bam" -b -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.bam"
            picard SortSam  {params.picard_Xmx}  INPUT="{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.bam" OUTPUT="{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=YES
            samtools view -h "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.bam" | grep -v -Pe "SN:chr[1234567890XY]+-" | samtools view -@1 -b -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
            samtools index -@1 "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
        fi

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



rule B15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part06__really_call_variants:
    input:
        lambda wildcards: ["result/S15_1__get_sample_RNA_editing_sites_v3/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/" + str(temp_row.INDEXER_PARAMETERS) + "/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part06__really_call_variants" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_B15_1__step02__part06=touch("result/B15_1__get_sample_RNA_editing_sites_v3/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part06__really_call_variants")


        
rule S15_1__get_sample_RNA_editing_sites_v3__step03__extract_bam_subset:
    input:
        flag_S15_1__step02__part05="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part05__locally_realign_reads"
    output:
        flag_S15_1__step03=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/recal.bam.subset/{REGION_ID}/finished.step03__extract_bam_subset"),
        recal_subset_bam_filename="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/recal.bam.subset/{REGION_ID}/recal.subset.bam"
    params:
        calling_result_directory="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/"
    threads:
        config['threads_subset_bam']
    shell:
        """
        region_id="{wildcards.REGION_ID}"

        if [[ $region_id == "BLCAP.Y2C.u10.d10" ]] ; then
            samtools view -@{threads} -b -o "{output.recal_subset_bam_filename}" "{params.calling_result_directory}/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.recal.bam" chr20:37519160-37519180
            samtools index -@{threads} "{output.recal_subset_bam_filename}"
        else
            echo "Unsupported REGION_ID: {wildcards.REGION_ID}" && exit 1;
        fi
        """



rule B15_1__get_sample_RNA_editing_sites_v3__step03__extract_bam_subset:
    input:
        lambda wildcards: ["result/S15_1__get_sample_RNA_editing_sites_v3/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/" + str(temp_row.INDEXER_PARAMETERS) + "/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/recal.bam.subset/{REGION_ID}/finished.step03__extract_bam_subset" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_B15_1__step03=touch("result/B15_1__get_sample_RNA_editing_sites_v3/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/recal.bam.subset/{REGION_ID}/finished.step03__extract_bam_subset")



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



rule S15_1__get_sample_RNA_editing_sites_v3____step07__apply_complex_filter____part09__reformat_data_as_standard_rich_vcf:
    input:
        flag_S15_1__step07__part07="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part07__separate_RepNOTAlu_and_nonRep"
    output:
        flag_S15_1__step07__part09=touch("result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/finished.step07__apply_complex_filter_1____part09__reformat_data_as_standard_rich_vcf")
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
            echo "##fileformat=VCFv4.2" > "{params.complex_filter_directory}/alignment.$subset".real.rich.vcf
            echo '##INFO=<ID=SUBSET,Number=1,Type=String,Description="Type of subset">'  >> "{params.complex_filter_directory}/alignment.$subset".real.rich.vcf
            echo '##FORMAT=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">' >> "{params.complex_filter_directory}/alignment.$subset".real.rich.vcf
            echo '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">'  >> "{params.complex_filter_directory}/alignment.$subset".real.rich.vcf
            echo '##FORMAT=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">'  >> "{params.complex_filter_directory}/alignment.$subset".real.rich.vcf
            echo "#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,{wildcards.SAMPLE}" | tr ',' '\t'  >> "{params.complex_filter_directory}/alignment.$subset".real.rich.vcf
            cat "{params.complex_filter_directory}/alignment.$subset".vcf | awk -v vSUBSET=$subset 'BEGIN{{OFS="\t"}}{{split($3, AN_AND_AC, ","); print $1,$2,$1"_"$2"_"$4"_"$5,$4,$5,".",".","SUBSET="vSUBSET,"AC:AN:AF",AN_AND_AC[2]":"AN_AND_AC[1]":"$6}}' >> "{params.complex_filter_directory}/alignment.$subset".real.rich.vcf
            rm -fr "{params.complex_filter_directory}/alignment.$subset".real.rich.vcf.gz*
            bgzip "{params.complex_filter_directory}/alignment.$subset".real.rich.vcf
            tabix "{params.complex_filter_directory}/alignment.$subset".real.rich.vcf.gz
        done

        bcftools concat --allow-overlaps -o "{params.complex_filter_directory}/alignment.all".real.rich.vcf -Ov "{params.complex_filter_directory}/alignment.Alu".real.rich.vcf.gz "{params.complex_filter_directory}/alignment.RepNOTAlu".real.rich.vcf.gz "{params.complex_filter_directory}/alignment.nonRep".real.rich.vcf.gz
        rm -fr "{params.complex_filter_directory}/alignment.all".real.rich.vcf.gz*
        bgzip -c "{params.complex_filter_directory}/alignment.all".real.rich.vcf > "{params.complex_filter_directory}/alignment.all".real.rich.vcf.gz
        bcftools index "{params.complex_filter_directory}/alignment.all".real.rich.vcf.gz
        """


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
        temp_vcf_path_collection = lambda wildcards: ["result/S15_1__get_sample_RNA_editing_sites_v3/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/"  + str(temp_row.INDEXER_PARAMETERS) + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/alignment.all.real.rich.vcf.gz" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)],
        result_directory="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{COMPLEX_FILTER}/{COMPLEX_FILTER_PARAMETERS}/"
    threads:
        config['threads_merging_vcfs']
    shell:
        """
        rm -fr {params.result_directory}/*
        mkdir -p {params.result_directory}
        echo {params.temp_vcf_path_collection} | tr ' ' '\n' > {params.result_directory}/vcf_path_to_merge.txt
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





rule S06_1__get_expression_level__step01__align_sample:
    input:
        flag_S01_2__step01="result/S01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/finished.step01__trim_RNA_Seq",
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


rule S06_1__get_expression_level__step02__call_expression:
    input:
        flag_S06_1__step01_2="result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01_2__merge_bams_of_runs_for_the_same_sample",
        alignment_merged_filename="result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/alignment.merged.bam"
    output:
        flag_S06_1__step02=touch("result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/finished.step02__call_expression")
    params:
        merge_result_directory="result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/",
        calling_expression_result_directory="result/S06_1__get_expression_level/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/",
        reference_GTF_filename="external/reference.gene.annotation/GENCODE.annotation/{GENCODE_VERSION}/gencode.annotation.gtf"
    threads:
        config['threads_calling_expression']
    shell:
        """
        rm -fr {params.calling_expression_result_directory}/*
        mkdir -p {params.calling_expression_result_directory}
        if [[ {wildcards.EXPRESSION_CALLER} == 'stringtie' ]]; then
            stringtie_parameters={wildcards.EXPRESSION_CALLER_PARAMETERS}
            if [[ $stringtie_parameters == 'none' ]]; then
                stringtie_parameters=""
            fi
            stringtie $stringtie_parameters -e -B -p {threads} -G {params.reference_GTF_filename} -o {params.calling_expression_result_directory}/assembly.gtf {input.alignment_merged_filename}
        else
            echo 'Unsupported EXPRESSION CALLER {wildcards.EXPRESSION_CALLER}' && exit 1
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
        "./scripts/S21_1__merge_phenotype_tables/run.r"

