import pandas
import re

for temp_key in ['threads_indexing', 'threads_trimming', 'threads_aligning', 'threads_calling_variants', 'threads_merging_bams', 'threads_re_faidx']:
    if temp_key not in config.keys():
        config[temp_key] = 20

for temp_key in ['threads_auxiliary_processing', 'threads_calling_variants_with_samtools_mpileup']:
    if temp_key not in config.keys():
        config[temp_key] = 1

if 'default_Xmx' not in config.keys():
    config['default_Xmx']='-Xmx60G'

for temp_key in ['picard_Xmx', 'GATK_BaseRecalibrator_Xmx', 'GATK_PrintReads_Xmx']:
    if temp_key not in config.keys():
        config[temp_key] = config['default_Xmx']


def __get_raw_fastq_path_by_TYPE_and_DATASET_and_SAMPLE_and_RUN(wildcards):
    parent_path = "external/RNA-Seq-with-Run/{DATASET}/{SAMPLE}/{RUN}/RNA/"
    type_end = re.sub("^([^-]+)-.*", "\\1", wildcards.TYPE)
    if type_end == "single":
        return parent_path + "r.fastq.gz"
    elif type_end == "paired":
        return [parent_path + "r1.fastq.gz", parent_path + "r2.fastq.gz"]

rule PS01_2__trim_RNA_Seq____step01__trim_RNA_Seq_by_fastp:
    input:
        input_read_filenames=__get_raw_fastq_path_by_TYPE_and_DATASET_and_SAMPLE_and_RUN
    output:
        flag_PS01_2__step01=touch("result/PS01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/finished.step01__trim_RNA_Seq_by_fastp")
    params:
        raw_fastq_parent_path="external/RNA-Seq-with-Run/{DATASET}/{SAMPLE}/{RUN}/RNA/",
        trim_result_directory="result/PS01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/"
    threads:
        config['threads_trimming']
    shell:
        """
        type_end=`echo {wildcards.TYPE} | cut -f 1 -d '-'`
        if [[ $type_end == 'paired' ]]; then
            if [[ "{wildcards.FASTP_PARAMETERS}" == 'filter-by-soapnuke' ]]; then
                tools/SOAPnuke/SOAPnuke filter -1 "{params.raw_fastq_parent_path}/r1.fastq.gz" -2 "{params.raw_fastq_parent_path}/r2.fastq.gz" -C "r1.fastp.fastq.gz" -D "r2.fastp.fastq.gz" -o "{params.trim_result_directory}" -T {threads}
            else
                echo "Unsupported FASTP_PARAMETERS: {wildcards.FASTP_PARAMETERS}" && exit 1
            fi
        elif [[ $type_end == 'single' ]]; then
            if [[ "{wildcards.FASTP_PARAMETERS}" == 'trim-15bp-off-3prime-and-filter-by-soapnuke' ]]; then
                fastp --disable_adapter_trimming --disable_trim_poly_g --disable_quality_filtering --disable_length_filtering --trim_tail1 15 -i {params.raw_fastq_parent_path}/r.fastq.gz -o "{params.trim_result_directory}/r.15bptrimmed.fastq.gz" -w {threads}
                tools/SOAPnuke/SOAPnuke filter -1 "{params.trim_result_directory}/r.15bptrimmed.fastq.gz" -C "r.fastp.fastq.gz" -o "{params.trim_result_directory}" -T {threads}
            else
                echo "Unsupported FASTP_PARAMETERS: {wildcards.FASTP_PARAMETERS}" && exit 1
            fi
        else
            echo "Unsupported type: {wildcards.TYPE}" && exit 1
        fi
        """
        
    
## in these series, "TYPE" is of the form "single-50" or "paired-75-60". The numbers are the avgSpotLen (read length) reported by NCBI or some other tools. Currently we ask the users to supply it instead of infer it from fastq files because it is cumbersome in the engineering sense. However, only the first part (single or paired) is used currently. The read length information is just for log purpose only.

## when "{wildcards.INDEXER}" == 'bwa-index-10.1038_nmeth.2330', the "{wildcards.INDEXER_PARAMETER}" is the read length
rule Ps05_1__index_contig_with_annotation:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}"
    output:
        flag_Ps05_1=touch("result/Ps05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/finished")
    params:
        reference_GTF_filename="external/reference.gene.annotation/GENCODE.annotation/{GENCODE_VERSION}/gencode.annotation.gtf",
        reference_UCSC_Table_Browser_knownGene_filename="./external/UCSC.Table.Browser.knownGene.GENCODE/{GENCODE_VERSION}/knownGene",
        contigs_index_directory="result/Ps05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/",
        contigs_index_prefix="result/Ps05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/prefix"
    threads:
        config['threads_indexing']
    conda:
        "py2.7.yaml"
    shell:
        """
        rm -fr "{params.contigs_index_directory}"/*
        mkdir -p "{params.contigs_index_directory}"
        if [[ "{wildcards.INDEXER}" == 'tophat2-index' ]]; then
            if [[ "{wildcards.INDEXER_PARAMETERS}" == 'tophat2-index-default' ]]; then
                bowtie2-build {input.contigs_fasta_filename} {params.contigs_index_directory}/genome-index
                ln -s -r {input.contigs_fasta_filename} {params.contigs_index_directory}/genome-index.fa
                if [[ "{wildcards.GENCODE_VERSION}" == "Ensembl75.GRCh37" ]] ; then
                    cat {params.reference_GTF_filename} | awk '{{ if ($1 ~ /^#/) {{print}} else {{print "chr"$0 }} }}' > {params.contigs_index_directory}/GTF.chr.corrected.gtf
                else 
                   ln -s -r {params.reference_GTF_filename} {params.contigs_index_directory}/GTF.chr.corrected.gtf
                fi
                tophat2 -G {params.contigs_index_directory}/GTF.chr.corrected.gtf --transcriptome-index={params.contigs_index_directory}/transcriptome-index --output-dir {params.contigs_index_directory}/output-dir {params.contigs_index_directory}/genome-index
            else
                echo "Indexer parameters {wildcards.INDEXER_PARAMETERS} not supported." && exit 1;
            fi
        else
            echo "Indexer {wildcards.INDEXER} not supported." && exit 1;
        fi
        """


rule PS15_1__get_sample_RNA_editing_sites_v3__step01__align_sample:
    input:
        flag_PS01_2__step01="result/PS01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/finished.step01__trim_RNA_Seq_by_fastp",
        flag_Ps05_1="result/Ps05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/finished"
    output:
        flag_PS15_1__step01=touch("result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01__align_sample")
    params:
        contigs_index_directory="result/Ps05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/",
        contigs_index_prefix="result/Ps05_1__index_contig_with_annotation/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/prefix",
        trim_result_directory="result/PS01_2__trim_RNA_Seq/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/",
        alignment_result_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/{RUN}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/"
    threads:
        config['threads_aligning']
    conda:
        "py2.7.yaml"
    shell:
        """
        type_end=`echo {wildcards.TYPE} | cut -f 1 -d '-'`
        rm -fr "{params.alignment_result_directory}/"*
        mkdir -p "{params.alignment_result_directory}/"
        if [[ "{wildcards.ALIGNER}" == 'tophat2' ]]; then
            aligner_parameters="{wildcards.ALIGNER_PARAMETERS}"
            if [[ $aligner_parameters != 'tophat2-default' ]]; then
                echo 'Unsupported ALIGNER_PARAMETERS {wildcards.ALIGNER_PARAMETERS}' && exit 1
            fi
            if [[ $type_end == single ]]; then
                tophat2 --read-gap-length 3 --read-edit-dist 3 --no-novel-juncs -no-novel-indels --num-threads {threads} --output-dir "{params.alignment_result_directory}/" --transcriptome-index={params.contigs_index_directory}/transcriptome-index {params.contigs_index_directory}/genome-index "{params.trim_result_directory}/r.fastp.fastq.gz"
                ln -s -r "{params.alignment_result_directory}/"accepted_hits.bam "{params.alignment_result_directory}/"alignment.bam
            elif [[ $type_end == paired ]]; then
                tophat2 --read-gap-length 3 --read-edit-dist 3 --no-novel-juncs -no-novel-indels --num-threads {threads} --output-dir "{params.alignment_result_directory}/" --transcriptome-index={params.contigs_index_directory}/transcriptome-index {params.contigs_index_directory}/genome-index "{params.trim_result_directory}/r1.fastp.fastq.gz" "{params.trim_result_directory}/r2.fastp.fastq.gz"
                ln -s -r "{params.alignment_result_directory}/"accepted_hits.bam "{params.alignment_result_directory}/"alignment.bam
            else
                echo 'Unsupported read type {wildcards.TYPE}' && exit 1
            fi
        else
            echo 'Unsupported ALIGNER {wildcards.ALIGNER}' && exit 1
        fi
        """



rule Ps05_3__index_contig_with_samtools_faidx:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}"
    output:
        flag_Ps05_3=touch("result/Ps05_3__index_contig_with_samtools_faidx/{CONTIGS_FASTA_FILENAME}/finished")
    threads:
        1
    shell:
        """
        samtools faidx {input.contigs_fasta_filename}
        """

rule Ps05_4__create_sequence_dictionary_for_contig_with_picard:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}"
    output:
        flag_Ps05_4=touch("result/Ps05_4__create_sequence_dictionary_for_contig_with_picard/{CONTIGS_FASTA_FILENAME}/finished")
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


rule PS15_1__get_sample_RNA_editing_sites_v3__step01_2__merge_bams_of_runs_for_the_same_sample:
    input:
        flags_collection_PS15_1__step01 = lambda wildcards: ["result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/" + temp_run + "/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01__align_sample" for temp_run in __get_RUNs_by_DATASET_and_SAMPLE(wildcards.DATASET, wildcards.SAMPLE)]
    output:
        flag_PS15_1__step01_2=touch("result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01_2__merge_bams_of_runs_for_the_same_sample"),
        alignment_merged_filename="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/alignment.merged.bam"
    params:
        alignment_bams_prefix = "result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/",
        alignment_bams_suffix = "/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/alignment.bam",
        alignment_bams_runs = lambda wildcards: __get_RUNs_by_DATASET_and_SAMPLE(wildcards.DATASET, wildcards.SAMPLE),
        number_of_alignments_to_merge = lambda wildcards: len(__get_RUNs_by_DATASET_and_SAMPLE(wildcards.DATASET, wildcards.SAMPLE)),
        merge_result_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/",
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



rule PS15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part01__sort_bam:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}",
        flag_Ps05_3="result/Ps05_3__index_contig_with_samtools_faidx/{CONTIGS_FASTA_FILENAME}/finished",
        flag_Ps05_4="result/Ps05_4__create_sequence_dictionary_for_contig_with_picard/{CONTIGS_FASTA_FILENAME}/finished",
        flag_PS15_1__step01_2="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/finished.step01_2__merge_bams_of_runs_for_the_same_sample"
    output:
        flag_PS15_1__step02__part01=touch("result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part01__sort_bam")
    params:
        alignment_merged_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/",
        calling_result_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        picard_Xmx=config["picard_Xmx"]
    threads:
        config['threads_calling_variants']
    shell:
        """
        mkdir -p "{params.calling_result_directory}/"
        picard SortSam {params.picard_Xmx} INPUT="{params.alignment_merged_directory}/alignment.merged.bam" OUTPUT="{params.calling_result_directory}/alignment.sorted.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
        """

rule PS15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part02__label_bam:
    input:
        flag_PS15_1__step02__part01="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part01__sort_bam"
    output:
        flag_PS15_1__step02__part02=touch("result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part02__label_bam")
    params:
        calling_result_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        picard_Xmx=config["picard_Xmx"]
    threads:
        config['threads_calling_variants']
    shell:
        """
        picard AddOrReplaceReadGroups  {params.picard_Xmx}  I="{params.calling_result_directory}/alignment.sorted.bam" O="{params.calling_result_directory}/alignment.sorted.withRG.bam" SO=coordinate RGID="{wildcards.SAMPLE}" RGLB="{wildcards.SAMPLE}" RGPL=illumina RGPU=illumina RGSM="{wildcards.SAMPLE}" VALIDATION_STRINGENCY=SILENT
        """


rule PS15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part03__mark_duplicates:
    input:
        flag_PS15_1__step02__part02="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part02__label_bam"
    output:
        flag_PS15_1__step02__part03=touch("result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part03__mark_duplicates")
    params:
        calling_result_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        picard_Xmx=config["picard_Xmx"]
    threads:
        config['threads_calling_variants']
    shell:
        """
        picard MarkDuplicates  {params.picard_Xmx}  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=2000 INPUT="{params.calling_result_directory}/alignment.sorted.withRG.bam" OUTPUT="{params.calling_result_directory}/alignment.sorted.withRG.dedup.bam" METRICS_FILE="{params.calling_result_directory}/alignment.sorted.withRG.dedup.matrix" CREATE_INDEX=true, VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
        """

rule Ps05_5__create_contig_info_in_karyotypic_order_for_picard_ReorderSam:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}",
        flag_Ps05_3="result/Ps05_3__index_contig_with_samtools_faidx/{CONTIGS_FASTA_FILENAME}/finished"
    output:
        flag_Ps05_4=touch("result/Ps05_5__create_contig_info_in_karyotypic_order_for_picard_ReorderSam/{CONTIGS_FASTA_FILENAME}/finished"),
        reordered_contig_fa_filename="result/Ps05_5__create_contig_info_in_karyotypic_order_for_picard_ReorderSam/{CONTIGS_FASTA_FILENAME}/reordered.contig.fa",
        reordered_contig_dict_filename="result/Ps05_5__create_contig_info_in_karyotypic_order_for_picard_ReorderSam/{CONTIGS_FASTA_FILENAME}/reordered.contig.dict"
    params:
        result_directory="result/Ps05_5__create_contig_info_in_karyotypic_order_for_picard_ReorderSam/{CONTIGS_FASTA_FILENAME}/"
    threads:
        config['threads_re_faidx']
    shell:
        """
        cut -f1 external/contigs/{wildcards.CONTIGS_FASTA_FILENAME}.fai | grep -v -P "(_|^chrM$)" | sort -k1V | parallel -j {threads} -k "samtools faidx {input.contigs_fasta_filename} {{}}" > {params.result_directory}/temp.contig.main.fa
        cut -f1 external/contigs/{wildcards.CONTIGS_FASTA_FILENAME}.fai | grep -P "(_|^chrM$)" | sort -k1V | parallel -j {threads} -k "samtools faidx {input.contigs_fasta_filename} {{}}" > {params.result_directory}/temp.contig.others.fa
        cat {params.result_directory}/temp.contig.main.fa {params.result_directory}/temp.contig.others.fa > {params.result_directory}/reordered.contig.fa
        samtools faidx {params.result_directory}/reordered.contig.fa
        picard CreateSequenceDictionary R={params.result_directory}/reordered.contig.fa O={params.result_directory}/reordered.contig.dict
        """


rule PS15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part04__reorder_sam:
    input:
        flag_PS15_1__step02__part03="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part03__mark_duplicates",
        flag_Ps05_4="result/Ps05_5__create_contig_info_in_karyotypic_order_for_picard_ReorderSam/{CONTIGS_FASTA_FILENAME}/finished",
        reordered_contig_dict_filename="result/Ps05_5__create_contig_info_in_karyotypic_order_for_picard_ReorderSam/{CONTIGS_FASTA_FILENAME}/reordered.contig.dict"
    output:
        flag_PS15_1__step02__part04=touch("result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part04__reorder_sam")
    params:
        calling_result_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        picard_Xmx=config["picard_Xmx"]
    threads:
        config['threads_calling_variants']
    shell:
        """
        picard ReorderSam  {params.picard_Xmx} -I "{params.calling_result_directory}/alignment.sorted.withRG.dedup.bam" -O "{params.calling_result_directory}/alignment.sorted.withRG.dedup.reordered.bam" -SD "{input.reordered_contig_dict_filename}"
        samtools index -@ {threads} "{params.calling_result_directory}/alignment.sorted.withRG.dedup.reordered.bam"
        """


rule Ps05_2__reformat_dbSNP_vcf:
    input:
        dbSNP_vcf_gz_filename="external/dbSNP.vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/dbSNP.vcf.gz"
    output:
        flag_Ps05_2=touch("result/Ps05_2__reformat_dbSNP_vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished"),
        dbSNP_vcf_reformatted_gz_filename="result/Ps05_2__reformat_dbSNP_vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/dbSNP.vcf.reformatted.gz"
    threads:
        config['threads_auxiliary_processing']
    shell:
        """
        zcat {input.dbSNP_vcf_gz_filename} | awk '{{if($1 ~ /^#/){{print $0}} else {{print "chr"$0}} }}' | bgzip -@{threads}  > {output.dbSNP_vcf_reformatted_gz_filename}
        tabix {output.dbSNP_vcf_reformatted_gz_filename}
        """

        
rule PS15_1__get_sample_RNA_editing_sites_v3__step02__call_variants____part05__recalibrate_base_quality:
    input:
        flag_PS15_1__step02__part03="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part04__reorder_sam",
        flag_Ps05_4="result/Ps05_5__create_contig_info_in_karyotypic_order_for_picard_ReorderSam/{CONTIGS_FASTA_FILENAME}/finished",
        reordered_contig_fa_filename="result/Ps05_5__create_contig_info_in_karyotypic_order_for_picard_ReorderSam/{CONTIGS_FASTA_FILENAME}/reordered.contig.fa",
        flag_Ps05_2="result/Ps05_2__reformat_dbSNP_vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished",
        dbSNP_vcf_reformatted_gz_filename="result/Ps05_2__reformat_dbSNP_vcf/{dbSNP_VERSION}/{dbSNP_SUBSET}/dbSNP.vcf.reformatted.gz"
    output:
        flag_PS15_1__step02__part05=touch("result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part05__recalibrate_base_quality")
    params:
        calling_result_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        GATK_BaseRecalibrator_Xmx=config["GATK_BaseRecalibrator_Xmx"],
        GATK_PrintReads_Xmx=config["GATK_PrintReads_Xmx"]
    threads:
        config['threads_calling_variants']
    conda:
        "jre.1.7.0.yaml"
    shell:
        """
        javaexec="java"
        if [[ "{wildcards.VARIANT_CALLER}" == "GATK-2.8-1" ]]; then
            javaexec="/gpfs2/gaog_pkuhpc/users/dingy/HERE-dev/.snakemake/conda/ed2ab993/x86_64-conda_cos6-linux-gnu/sysroot/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.131.x86_64/jre/bin/java"
        fi
        $javaexec {params.GATK_BaseRecalibrator_Xmx} -jar ./tools/{wildcards.VARIANT_CALLER}/GenomeAnalysisTK.jar -R "{input.reordered_contig_fa_filename}" -T BaseRecalibrator -I "{params.calling_result_directory}/alignment.sorted.withRG.dedup.reordered.bam" -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.recal_data.grp" -knownSites {input.dbSNP_vcf_reformatted_gz_filename} -U ALLOW_N_CIGAR_READS
        $javaexec {params.GATK_PrintReads_Xmx} -jar ./tools/{wildcards.VARIANT_CALLER}/GenomeAnalysisTK.jar -R "{input.reordered_contig_fa_filename}" -T PrintReads -I "{params.calling_result_directory}/alignment.sorted.withRG.dedup.reordered.bam" -o "{params.calling_result_directory}/alignment.sorted.withRG.dedup.reordered.recal.bam" -BQSR "{params.calling_result_directory}/alignment.sorted.withRG.dedup.recal_data.grp" -U ALLOW_N_CIGAR_READS
        """

