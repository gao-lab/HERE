import pandas
import os

for temp_key in ['threads_generate_vcf', 'threads_annotate_with_snpEff', 'threads_generate_bed', 'threads_pigz_TargetScan_output']:
    if temp_key not in config.keys():
        config[temp_key] = 20

for temp_key in ['snpEff_Xmx']:
    if temp_key not in config.keys():
        config[temp_key] = "-Xmx150G"



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







rule A02_4__check_fine_recurrence_profile_for_a_subset_of_samples:
    input:
        flag_S52_3="result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_merged_with_coverage_dt_txt_gz_filename="result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt.txt.gz",
        flag_S51_6="result/S51_6__get_snpEff_annotation_subset_of_filtered_result/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_variant_only_snpEff_annotation_dt_txt_gz_filename="result/S51_6__get_snpEff_annotation_subset_of_filtered_result/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt.txt.gz",
        total_sample_count_for_normal_stages_dt_csv_filename="report.ver2/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/total.sample.count.for.normal.stages.dt.csv",
        total_sample_count_for_GSE133854_all_dt_csv_filename="report.ver2/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/total.sample.count.for.GSE133854.all.dt.csv"
    output:
        flag_A02_4=touch("result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/finished"),
        subset_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.dt.txt.gz",
        subset_site_recurrence_comparison_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.site.recurrence.comparison.dt.txt.gz",
        subset_site_recurrence_comparison_recurrent_edits_only_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.site.recurrence.comparison.recurrent.edits.only.dt.txt.gz",
        subset_recurrent_edits_only_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.recurrent.edits.only.dt.txt.gz",
        subset_site_recurrence_comparison_CJ_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.site.recurrence.comparison.CJ.dt.txt.gz",
        snpEff_annotation_for_subset_recurrent_edits_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/snpEff.annotation.for.subset.recurrent.edits.dt.txt.gz",
        snpEff_annotation_for_subset_recurrent_edits_on_valid_transcripts_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt.txt.gz",
        snpEff_annotation_for_subset_recurrent_edits_on_valid_genes_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt.txt.gz",
        subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz",
        snpEff_annotation_for_subset_observed_edits_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/snpEff.annotation.for.subset.observed.edits.dt.txt.gz",
        snpEff_annotation_for_subset_observed_edits_on_valid_transcripts_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/snpEff.annotation.for.subset.observed.edits.on.valid.transcripts.dt.txt.gz",
        snpEff_annotation_for_subset_observed_edits_on_valid_genes_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/snpEff.annotation.for.subset.observed.edits.on.valid.genes.dt.txt.gz",
        subset_observed_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.observed.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz"
    script:
        "./scripts/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/run.R"


        

rule sA02_5__generate_all_transcribable_Alu_A_vcf____step01__generate_Alu_only_vcf:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}",
        repeatmasker_repFamily_Alu_bed_filename="external/UCSC.Table.Browser.repeatmasker/repFamily.Alu/repeatmasker.bed"
    output:
        flag_sA02_5_step01=touch("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/finished.step01__generate_Alu_only_vcf")
    params:
        result_directory="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/"
    threads:
        config['threads_generate_vcf']
    script:
        "./scripts/sA02_5__generate_all_transcribable_Alu_A_vcf/step01__generate_Alu_only_vcf/run.R"


rule sA02_5__generate_all_transcribable_Alu_A_vcf____step02__annotate_a_single_Alu_only_vcf_with_snpEff:
    input:
        flag_sA02_5_step01="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/finished.step01__generate_Alu_only_vcf"
    output:
        flag_sA02_5_step02=touch("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step02__annotate_a_single_Alu_only_vcf_with_snpEff.{CHROMOSOME}"),
        single_contig_Alu_only_snpEff_annotated_vcf_gz_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.vcf.gz"
    params:
        single_contig_Alu_only_vcf_gz_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{CHROMOSOME}.Alu.only.vcf.gz",
        snpEff_config_filename="result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/snpEff.config",
        snpEff_database_name="{CONTIGS_FASTA_FILENAME}.GENCODE.{GENCODE_VERSION}",
        snpEff_Xmx=config['snpEff_Xmx']
    threads:
        config['threads_annotate_with_snpEff']
    shell:
        """
        snpEff ann {params.snpEff_Xmx} -verbose -no-intergenic -no-upstream -no-downstream -config {params.snpEff_config_filename} {params.snpEff_database_name} {params.single_contig_Alu_only_vcf_gz_filename} | bcftools view -o {output.single_contig_Alu_only_snpEff_annotated_vcf_gz_filename} -Oz
        sleep 1
        tabix {output.single_contig_Alu_only_snpEff_annotated_vcf_gz_filename}
        """


rule sA02_5__generate_all_transcribable_Alu_A_vcf____step03__extract_annotated_records_from_a_single_Alu_only_vcf:
    input:
        flag_sA02_5_step02="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step02__annotate_a_single_Alu_only_vcf_with_snpEff.{CHROMOSOME}",
        single_contig_Alu_only_snpEff_annotated_vcf_gz_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.vcf.gz"
    output:
        flag_sA02_5_step03=touch("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step03__extract_annotated_records_from_a_single_Alu_only_vcf.{CHROMOSOME}"),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_vcf_gz_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.vcf.gz"
    threads:
        config['threads_generate_vcf']
    shell:
        """
        bcftools view --threads {threads} --output-file {output.single_contig_Alu_only_snpEff_annotated_annotated_records_only_vcf_gz_filename} -Oz --include 'INFO/ANN ~ "protein_coding.*A>G"' {input.single_contig_Alu_only_snpEff_annotated_vcf_gz_filename}
        sleep 1
        tabix {output.single_contig_Alu_only_snpEff_annotated_annotated_records_only_vcf_gz_filename}
        """


rule sA02_5__generate_all_transcribable_Alu_A_vcf____step04__split_a_single_annotated_records_of_Alu_only_vcf_into_subsets:
    input:
        flag_sA02_5_step03="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step03__extract_annotated_records_from_a_single_Alu_only_vcf.{CHROMOSOME}",
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_vcf_gz_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.vcf.gz"
    output:
        flag_sA02_5_step04=touch("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step04__split_a_single_annotated_records_of_Alu_only_vcf_into_subsets.{CHROMOSOME}"),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_intron_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.intron.bed",
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.3.prime.UTR.bed",
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.5.prime.UTR.bed",
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_CDS_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.CDS.bed"
    threads:
        config['threads_generate_bed']
    script:
        "./scripts/sA02_5__generate_all_transcribable_Alu_A_vcf/step04__split_a_single_annotated_records_of_Alu_only_vcf_into_subsets/run.R"


temp_chr_collection = ['chr' + str(temp_item) for temp_item in list(range(1, 23)) + ['X', 'Y']]

rule sA02_5__generate_all_transcribable_Alu_A_vcf____step05__combine_all_subsets:
    input:
        flag_sA02_5_step04_collection=expand("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step04__split_a_single_annotated_records_of_Alu_only_vcf_into_subsets.{CHROMOSOME}", CHROMOSOME=temp_chr_collection, allow_missing=True),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_intron_bed_filename_collection=expand("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.intron.bed", CHROMOSOME=temp_chr_collection, allow_missing=True),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_bed_filename_collection=expand("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.3.prime.UTR.bed", CHROMOSOME=temp_chr_collection, allow_missing=True),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_bed_filename_collection=expand("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.5.prime.UTR.bed", CHROMOSOME=temp_chr_collection, allow_missing=True),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_CDS_bed_filename_collection=expand("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.CDS.bed", CHROMOSOME=temp_chr_collection, allow_missing=True)
    output:
        flag_sA02_5_step04=touch("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step05__combine_all_subsets"),
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_intron_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.intron.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.3.prime.UTR.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.5.prime.UTR.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_CDS_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.CDS.combined.bed"
    threads:
        1
    shell:
        """
        cat {input.single_contig_Alu_only_snpEff_annotated_annotated_records_only_intron_bed_filename_collection} > {output.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_intron_combined_bed_filename}
        cat {input.single_contig_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_bed_filename_collection} > {output.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_combined_bed_filename}
        cat {input.single_contig_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_bed_filename_collection} > {output.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_combined_bed_filename}
        cat {input.single_contig_Alu_only_snpEff_annotated_annotated_records_only_CDS_bed_filename_collection} > {output.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_CDS_combined_bed_filename}
        """



rule A02_5__check_transcribable_Alu_A_coverage_of_merged_bam:
    input:
        flag_sA02_5_step04="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/finished.step05__combine_all_subsets",
        ## all_contigs_Alu_only_snpEff_annotated_annotated_records_only_intron_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.intron.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.3.prime.UTR.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.5.prime.UTR.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_CDS_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.CDS.combined.bed",
        ## flag_S15_1__step02__part04="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/{INDEXER_PARAMETERS}/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/finished.step02__call_variants____part04__correct_coordinates",
        alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/{INDEXER_PARAMETERS}/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
    output:
        flag_A02_5=touch("result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished"),
        ## alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.intron.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_5_prime_UTR_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.5.prime.UTR.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_3_prime_UTR_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.3.prime.UTR.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_CDS_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.CDS.bed.txt.gz"
    threads:
        1
    shell:
        """
        samtools depth -d 1 -b {input.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_combined_bed_filename} {input.alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename} | gzip > {output.alignment_merged_depth_on_transcribable_Alu_A_5_prime_UTR_bed_txt_gz_filename}
        samtools depth -d 1 -b {input.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_combined_bed_filename} {input.alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename} | gzip > {output.alignment_merged_depth_on_transcribable_Alu_A_3_prime_UTR_bed_txt_gz_filename}
        samtools depth -d 1 -b {input.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_CDS_combined_bed_filename} {input.alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename} | gzip > {output.alignment_merged_depth_on_transcribable_Alu_A_CDS_bed_txt_gz_filename}
        """



rule BA02_5__check_transcribable_Alu_A_coverage_of_merged_bam:
    input:
        lambda wildcards: ["result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/finished" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_BA02_5=touch("result/BA02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished")




rule A02_5__check_transcribable_Alu_A_coverage_of_merged_bam____patch01__intron:
    input:
        flag_sA02_5_step04="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/finished.step05__combine_all_subsets",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_intron_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.intron.combined.bed",
        ## flag_S15_1__step02__part04="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/{INDEXER_PARAMETERS}/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/finished.step02__call_variants____part04__correct_coordinates",
        alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/{INDEXER_PARAMETERS}/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
    output:
        flag_A02_5_patch01=touch("result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished.patch01__intron"),
        alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.intron.bed.txt.gz"
    threads:
        1
    shell:
        """
        samtools depth -d 1 -b {input.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_intron_combined_bed_filename} {input.alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename} | gzip > {output.alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename}
        """

rule BA02_5__check_transcribable_Alu_A_coverage_of_merged_bam____patch01__intron:
    input:
        lambda wildcards: ["result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/finished.patch01__intron" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_B52_1=touch("result/BA02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished.patch01__intron")


rule A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam:
    input:
        flag_A02_5="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished",
        ## alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.intron.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_5_prime_UTR_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.5.prime.UTR.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_3_prime_UTR_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.3.prime.UTR.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_CDS_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.CDS.bed.txt.gz"
    output:
        flag_A02_6=touch("result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished"),
        count_5_prime_UTR_txt_filename="result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/count.5.prime.UTR.txt",
        count_3_prime_UTR_txt_filename="result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/count.3.prime.UTR.txt",
        count_CDS_txt_filename="result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/count.CDS.txt"
    threads:
        1
    shell:
        """
        zcat {input.alignment_merged_depth_on_transcribable_Alu_A_5_prime_UTR_bed_txt_gz_filename} | awk '{{if ($3 > 0) print}}' | wc -l > {output.count_5_prime_UTR_txt_filename}
        zcat {input.alignment_merged_depth_on_transcribable_Alu_A_3_prime_UTR_bed_txt_gz_filename} | awk '{{if ($3 > 0) print}}' | wc -l > {output.count_3_prime_UTR_txt_filename}
        zcat {input.alignment_merged_depth_on_transcribable_Alu_A_CDS_bed_txt_gz_filename} | awk '{{if ($3 > 0) print}}' | wc -l > {output.count_CDS_txt_filename}
        """

rule BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam:
    input:
        flags_collection=lambda wildcards: ["result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/finished" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_BA02_5=touch("result/BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished")



rule A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam____patch01__intron:
    input:
        flag_A02_5="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished",
        alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.intron.bed.txt.gz"
    output:
        flag_A02_6__patch01=touch("result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished.patch01__intron"),
        count_intron_txt_filename="result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/count.intron.txt"
    threads:
        1
    shell:
        """
        zcat {input.alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename} | awk '{{if ($3 > 0) print}}' | wc -l > {output.count_intron_txt_filename}
        """

rule BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam____patch01__intron:
    input:
        flags_collection=lambda wildcards: ["result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/finished.patch01__intron" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_BA02_6__patch01=touch("result/BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished.patch01__intron")


rule A02_7__combine_summaries_of_distribution_of_transcribable_Alu_A_coverage_of_merged_bam:
    input:
        flag_BA02_6="result/BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        flag_BA02_6__patch01="result/BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished.patch01__intron"
    output:
        flag_A02_7=touch("result/A02_7__combine_summaries_of_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished"),
        final_count_csv_filename="result/A02_7__combine_summaries_of_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/final.count.csv"
    params:
        lambda wildcards: [ [temp_row.DATASET_NAME, temp_row.SAMPLE_NAME] + ["result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/" + wildcards.DATASET_COLLECTION_NAME + "/" + wildcards.DATASET_PHENOTYPE_COLLECTION_NAME + "/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/count." + Annotation_subset + ".txt" for Annotation_subset in ["5.prime.UTR", "3.prime.UTR", "CDS", "intron"]] for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    run:
        dataset_and_sample_and_counts_collection=params[0]
        all_records_list = []
        used_samples = []
        for dataset, sample, count_5_prime_UTR_filename, count_3_prime_UTR_filename, count_CDS_filename, count_intron_filename in dataset_and_sample_and_counts_collection:
            if sample in used_samples:
                print("The sample " + sample + " has been used")
                continue
            used_samples.append(sample)
            temp_counts_list = []
            for temp_filename in [count_5_prime_UTR_filename, count_3_prime_UTR_filename, count_CDS_filename, count_intron_filename]:
                temp_count = "NA"
                with open(temp_filename, "r") as temp_f:
                    temp_count=temp_f.readline().strip()
                temp_counts_list.append(temp_count)
            all_records_list.append([dataset, sample] + temp_counts_list)
        final_pd = pandas.DataFrame(all_records_list)
        final_pd.columns = ["dataset", "name", "count.5.prime.UTR", "count.3.prime.UTR", "count.CDS", "count.intron"]
        final_pd.to_csv(output.final_count_csv_filename)


rule A02_8__get_editing_effect_on_miRNA_binding_sites____step01__get_3UTR_blocks:
    input:
        reference_GTF_filename="external/reference.gene.annotation/GENCODE.annotation/{GENCODE_VERSION}/gencode.annotation.gtf"
    output:
        A02_8__step01_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step01__get_3UTR_blocks/{GENCODE_VERSION}/finished"),
        gencode_3utr_dt_csv_gz_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step01__get_3UTR_blocks/{GENCODE_VERSION}/gencode.3utr.dt.csv.gz",
        gencode_3utr_gff_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step01__get_3UTR_blocks/{GENCODE_VERSION}/gencode.3utr.gff"
    script:
        "./scripts/A02_8__get_editing_effect_on_miRNA_binding_sites/step01__get_3UTR_blocks/run.R"

rule A02_8__get_editing_effect_on_miRNA_binding_sites____step02__get_3UTR_bed_per_chromosome:
    input:
        gencode_3utr_gff_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step01__get_3UTR_blocks/{GENCODE_VERSION}/gencode.3utr.gff"
    output:
        A02_8__step02_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step02__get_3UTR_bed_per_chromosome/{GENCODE_VERSION}/finished"),
        gencode_3utr_bed_corrected_all_single_chromosome_filenames_list = ["result/A02_8__get_editing_effect_on_miRNA_binding_sites/step02__get_3UTR_bed_per_chromosome/{GENCODE_VERSION}/step2.gencode.3utr." + "chr" + str(chr) + ".col7.corrected.bed" for chr in list(range(1, 23)) + ["X", "Y"]]
    params:
        gencode_3utr_genePred_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step02__get_3UTR_bed_per_chromosome/{GENCODE_VERSION}/gencode.3utr.genePred",
        gencode_3utr_bed_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step02__get_3UTR_bed_per_chromosome/{GENCODE_VERSION}/gencode.3utr.bed",
        result_directory="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step02__get_3UTR_bed_per_chromosome/{GENCODE_VERSION}/"
    shell:
        """
        gtfToGenePred {input.gencode_3utr_gff_filename} {params.gencode_3utr_genePred_filename}
        genePredToBed {params.gencode_3utr_genePred_filename} {params.gencode_3utr_bed_filename}
        cat {params.gencode_3utr_bed_filename} | cut -f 1 | uniq | while read chr
        do 
            grep -P "^$chr\t" {params.gencode_3utr_bed_filename} | awk -v OFS='\t' '{{$7=$2; print}}' > {params.result_directory}/step2.gencode.3utr.$chr.col7.corrected.bed
        done
        """



rule A02_8__get_editing_effect_on_miRNA_binding_sites____step03__get_maf_fasta_per_chromosome:
    input:
        gencode_3utr_bed_corrected_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step02__get_3UTR_bed_per_chromosome/{GENCODE_VERSION}/step2.gencode.3utr.{CHR}.col7.corrected.bed",
        maf_for_a_single_chromosome_filename="external/UCSC.maf30way/{CHR}.maf",
        maf_index_for_a_single_chromosome_filename="external/UCSC.maf30way/{CHR}.maf.index"
    output:
        A02_8__step03_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step03__get_maf_fasta_per_chromosome/{GENCODE_VERSION}/finished.{CHR}"),
        gencode_3utr_maf_fasta_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step03__get_maf_fasta_per_chromosome/{GENCODE_VERSION}/gencode.3utr.{CHR}.maf.fasta"
    shell:
        """
        python ./tools/galaxy/interval_maf_to_merged_fasta.py -d hg38 -i {input.gencode_3utr_bed_corrected_for_a_single_chromosome_filename}  -m {input.maf_for_a_single_chromosome_filename} -o {output.gencode_3utr_maf_fasta_for_a_single_chromosome_filename} -G -I {input.maf_index_for_a_single_chromosome_filename}  --overwrite_with_gaps=TRUE --mafSourceType=user
        """


rule A02_8__get_editing_effect_on_miRNA_binding_sites____step04__generate_TargetScan_input_per_chromosome:
    input:
        gencode_3utr_maf_fasta_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step03__get_maf_fasta_per_chromosome/{GENCODE_VERSION}/gencode.3utr.{CHR}.maf.fasta",
        hg38_30way_taxonomy_filename="external/UCSC.maf30way/hg38.30way.taxonomy"
    output:
        A02_8__step04_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step04__generate_TargetScan_input_per_chromosome/{GENCODE_VERSION}/finished.{CHR}"),
        gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step04__generate_TargetScan_input_per_chromosome/{GENCODE_VERSION}/gencode.3utr.{CHR}.maf.table.with.species.id.txt"
    params:
        gencode_3utr_maf_table_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step04__generate_TargetScan_input_per_chromosome/{GENCODE_VERSION}/gencode.3utr.{CHR}.maf.table"
    shell:
        """
        cat {input.gencode_3utr_maf_fasta_for_a_single_chromosome_filename} | grep -v "^[ ]*$" | paste - - | sed 's/^>//' | sed 's/\./\t/' > {params.gencode_3utr_maf_table_for_a_single_chromosome_filename}
        Rscript ./scripts/A02_8__get_editing_effect_on_miRNA_binding_sites/step04__generate_TargetScan_input_per_chromosome/add.species.id.R {params.gencode_3utr_maf_table_for_a_single_chromosome_filename} {input.hg38_30way_taxonomy_filename} {output.gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename}
        """



rule A02_8__get_editing_effect_on_miRNA_binding_sites____step05__run_TargetScan_per_chromosome:
    input:
        gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step04__generate_TargetScan_input_per_chromosome/{GENCODE_VERSION}/gencode.3utr.{CHR}.maf.table.with.species.id.txt",
        miRNA_Family_Info_human_txt_filename="external/TargetScan/miR_Family_Info.human.txt"
    output:
        A02_8__step05_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step05__run_TargetScan_per_chromosome/{GENCODE_VERSION}/finished.{CHR}"),
        TargetScan_output_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step05__run_TargetScan_per_chromosome/{GENCODE_VERSION}/gencode.3utr.{CHR}.TargetScan.output"
    shell:
        """
        rm -f {output.TargetScan_output_filename}
        perl ./tools/targetscan_70/targetscan_70.pl {input.miRNA_Family_Info_human_txt_filename} {input.gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename} {output.TargetScan_output_filename}
        """


rule A02_8__get_editing_effect_on_miRNA_binding_sites____step06__compute_edit_relative_position_on_3UTR:
    input:
        gencode_3utr_dt_csv_gz_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step01__get_3UTR_blocks/{GENCODE_VERSION}/gencode.3utr.dt.csv.gz",
        snpEff_annotation_for_subset_recurrent_edits_on_valid_transcripts_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt.txt.gz"
    output:
        A02_8__step06_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step06__compute_edit_relative_position_on_3UTR/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished"),
        gencode_3utr_and_edit_CJ_dt_csv_gz_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step06__compute_edit_relative_position_on_3UTR/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/gencode.3utr.and.edit.CJ.dt.csv.gz"
    script:
        "./scripts/A02_8__get_editing_effect_on_miRNA_binding_sites/step06__compute_edit_relative_position_on_3UTR/run.R"


        

rule A02_8__get_editing_effect_on_miRNA_binding_sites____step07__get_edited_TargetScan_input_per_chromosome:
    input:
        gencode_3utr_and_edit_CJ_dt_csv_gz_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step06__compute_edit_relative_position_on_3UTR/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/gencode.3utr.and.edit.CJ.dt.csv.gz",
        gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step04__generate_TargetScan_input_per_chromosome/{GENCODE_VERSION}/gencode.3utr.{CHR}.maf.table.with.species.id.txt"
    output:
        A02_8__step07_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step07__get_edited_TargetScan_input_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished.{CHR}"),
        edited_gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step07__get_edited_TargetScan_input_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/edited.gencode.3utr.{CHR}.maf.table.with.species.id.txt"
    script:
        "./scripts/A02_8__get_editing_effect_on_miRNA_binding_sites/step07__get_edited_TargetScan_input_per_chromosome/run.R"


rule A02_8__get_editing_effect_on_miRNA_binding_sites____step08__run_TargetScan_on_edited_per_chromosome:
    input:
        edited_gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step07__get_edited_TargetScan_input_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/edited.gencode.3utr.{CHR}.maf.table.with.species.id.txt",
        miRNA_Family_Info_human_txt_filename="external/TargetScan/miR_Family_Info.human.txt"
    output:
        A02_8__step05_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step08__run_TargetScan_on_edited_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished.{CHR}"),
        TargetScan_on_edited_output_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step08__run_TargetScan_on_edited_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/edited.gencode.3utr.{CHR}.TargetScan.output"
    shell:
        """
        rm -f {output.TargetScan_on_edited_output_filename}
        perl ./tools/targetscan_70/targetscan_70.pl {input.miRNA_Family_Info_human_txt_filename} {input.edited_gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename} {output.TargetScan_on_edited_output_filename}
        """




rule A02_8__get_editing_effect_on_miRNA_binding_sites____step09__concatenate_TargetScan_results_across_all_chromosomes:
    input:
        TargetScan_output_filenames_list=["result/A02_8__get_editing_effect_on_miRNA_binding_sites/step05__run_TargetScan_per_chromosome/{GENCODE_VERSION}/gencode.3utr." + "chr" + str(chr) + ".TargetScan.output"  for chr in list(range(1, 23)) + ["X", "Y"]]
    output:
        A02_8__step09_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step09__concatenate_TargetScan_results_across_all_chromosomes/{GENCODE_VERSION}/finished"),
        concatenated_headless_TargetScan_output_gz_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step09__concatenate_TargetScan_results_across_all_chromosomes/{GENCODE_VERSION}/gencode.3utr.all.chromosomes.concatenated.headless.TargetScan.output.gz"
    threads:
        config['threads_pigz_TargetScan_output']
    shell:
        """
        for file in {input.TargetScan_output_filenames_list}
        do
            tail -n +2 $file
        done | pigz -c -p {threads} > {output.concatenated_headless_TargetScan_output_gz_filename}
        """

## no chrY
rule A02_8__get_editing_effect_on_miRNA_binding_sites____step10__concatenate_edited_TargetScan_results_across_all_chromosomes:
    input:
        edited_TargetScan_output_filenames_list=["result/A02_8__get_editing_effect_on_miRNA_binding_sites/step08__run_TargetScan_on_edited_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/edited.gencode.3utr." + "chr" + str(chr) + ".TargetScan.output"  for chr in list(range(1, 23)) + ["X"]]
    output:
        A02_8__step10_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step10__concatenate_edited_TargetScan_results_across_all_chromosomes/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished"),
        concatenated_headless_edited_TargetScan_output_gz_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step10__concatenate_edited_TargetScan_results_across_all_chromosomes/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/edited.gencode.3utr.all.chromosomes.but.chrY.concatenated.headless.TargetScan.output.gz"
    threads:
        config['threads_pigz_TargetScan_output']
    shell:
        """
        for file in {input.edited_TargetScan_output_filenames_list}
        do
            tail -n +2 $file
        done | pigz -c -p {threads} > {output.concatenated_headless_edited_TargetScan_output_gz_filename}
        """




rule A02_8__get_editing_effect_on_miRNA_binding_sites____step16__compute_all_edit_relative_position_on_3UTR:
    input:
        gencode_3utr_dt_csv_gz_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step01__get_3UTR_blocks/{GENCODE_VERSION}/gencode.3utr.dt.csv.gz",
        snpEff_annotation_for_subset_observed_edits_on_valid_transcripts_dt_txt_gz_filename="result/A02_4__check_fine_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/snpEff.annotation.for.subset.observed.edits.on.valid.transcripts.dt.txt.gz"
    output:
        A02_8__step16_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step16__compute_all_edit_relative_position_on_3UTR/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished"),
        gencode_3utr_and_all_edit_CJ_dt_csv_gz_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step16__compute_all_edit_relative_position_on_3UTR/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/gencode.3utr.and.all.edit.CJ.dt.csv.gz"
    script:
        "./scripts/A02_8__get_editing_effect_on_miRNA_binding_sites/step16__compute_all_edit_relative_position_on_3UTR/run.R"


        

rule A02_8__get_editing_effect_on_miRNA_binding_sites____step17__get_all_edited_TargetScan_input_per_chromosome:
    input:
        gencode_3utr_and_all_edit_CJ_dt_csv_gz_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step16__compute_all_edit_relative_position_on_3UTR/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/gencode.3utr.and.all.edit.CJ.dt.csv.gz",
        gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step04__generate_TargetScan_input_per_chromosome/{GENCODE_VERSION}/gencode.3utr.{CHR}.maf.table.with.species.id.txt"
    output:
        A02_8__step07_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step17__get_all_edited_TargetScan_input_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished.{CHR}"),
        all_edited_gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step17__get_all_edited_TargetScan_input_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/all.edited.gencode.3utr.{CHR}.maf.table.with.species.id.txt"
    script:
        "./scripts/A02_8__get_editing_effect_on_miRNA_binding_sites/step17__get_all_edited_TargetScan_input_per_chromosome/run.R"


rule A02_8__get_editing_effect_on_miRNA_binding_sites____step18__run_TargetScan_on_all_edited_per_chromosome:
    input:
        all_edited_gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step17__get_all_edited_TargetScan_input_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/all.edited.gencode.3utr.{CHR}.maf.table.with.species.id.txt",
        miRNA_Family_Info_human_txt_filename="external/TargetScan/miR_Family_Info.human.txt"
    output:
        A02_8__step18_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step18__run_TargetScan_on_all_edited_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished.{CHR}"),
        TargetScan_on_all_edited_output_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step18__run_TargetScan_on_all_edited_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/all.edited.gencode.3utr.{CHR}.TargetScan.output"
    shell:
        """
        rm -f {output.TargetScan_on_all_edited_output_filename}
        perl ./tools/targetscan_70/targetscan_70.pl {input.miRNA_Family_Info_human_txt_filename} {input.all_edited_gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename} {output.TargetScan_on_all_edited_output_filename}
        """


rule A02_8__get_editing_effect_on_miRNA_binding_sites____step20__concatenate_all_edited_TargetScan_results_across_all_chromosomes:
    input:
        all_edited_TargetScan_output_filenames_list=["result/A02_8__get_editing_effect_on_miRNA_binding_sites/step18__run_TargetScan_on_all_edited_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/all.edited.gencode.3utr." + "chr" + str(chr) + ".TargetScan.output"  for chr in list(range(1, 23)) + ["X", "Y"]]
    output:
        A02_8__step10_flag=touch("result/A02_8__get_editing_effect_on_miRNA_binding_sites/step20__concatenate_all_edited_TargetScan_results_across_all_chromosomes/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished"),
        concatenated_headless_all_edited_TargetScan_output_gz_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step20__concatenate_all_edited_TargetScan_results_across_all_chromosomes/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/all.edited.gencode.3utr.all.chromosomes.concatenated.headless.TargetScan.output.gz"
    threads:
        config['threads_pigz_TargetScan_output']
    shell:
        """
        for file in {input.all_edited_TargetScan_output_filenames_list}
        do
            tail -n +2 $file
        done | pigz -c -p {threads} > {output.concatenated_headless_all_edited_TargetScan_output_gz_filename}
        """



rule A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda____step00__get_miRNA_Family_human_fasta:
    input:
        miR_Family_Info_txt="external/TargetScan/miR_Family_Info.txt"
    output:
        A02_9__step00_flag=touch("result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step00__get_miRNA_Family_human_fasta/finished"),
        miR_Family_human_fasta="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step00__get_miRNA_Family_human_fasta/miR_Family_human.fasta"
    shell:
        """
        bash ./scripts/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step00__get_miRNA_Family_human_fasta/convert.TargetScan.Family.to.human.fasta.sh {input.miR_Family_Info_txt} > {output.miR_Family_human_fasta}
        """

        
rule A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda____step01__get_3UTR_sequences_from_maf_blocks:
    input:
        gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step04__generate_TargetScan_input_per_chromosome/{GENCODE_VERSION}/gencode.3utr.{CHR}.maf.table.with.species.id.txt"
    output:
        A02_9__step01_flag=touch("result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step01__get_3UTR_sequences_from_maf_blocks/{GENCODE_VERSION}/finished.{CHR}"),
        gencode_3utr_fasta_for_a_single_chromosome="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step01__get_3UTR_sequences_from_maf_blocks/{GENCODE_VERSION}/gencode.3utr.{CHR}.fasta"
    shell:
        """
        bash ./scripts/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step01__get_3UTR_sequences_from_maf_blocks/convert_maf_to_human_fasta.sh {input.gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename} > {output.gencode_3utr_fasta_for_a_single_chromosome}
        """

rule A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda____step02__run_miRanda_per_chromosome:
    input:
        gencode_3utr_fasta_for_a_single_chromosome="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step01__get_3UTR_sequences_from_maf_blocks/{GENCODE_VERSION}/gencode.3utr.{CHR}.fasta",
        miR_Family_human_fasta="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step00__get_miRNA_Family_human_fasta/miR_Family_human.fasta"
    output:
        A02_9__step02_flag=touch("result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step02__run_miRanda_per_chromosome/{GENCODE_VERSION}/finished.{CHR}"),
        gencode_3utr_miRanda_output_for_a_single_chromosome="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step02__run_miRanda_per_chromosome/{GENCODE_VERSION}/gencode.3utr.{CHR}.miRanda.output.gz"
    params:
        output_directory="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step02__run_miRanda_per_chromosome/{GENCODE_VERSION}/"
    shell:
        """
        mkdir -p {params.output_directory}
        rm -f {params.output_directory}/temp.{wildcards.CHR}.fasta
        ln -s -r {input.gencode_3utr_fasta_for_a_single_chromosome} {params.output_directory}/temp.{wildcards.CHR}.fasta
        ./tools/miRanda-1.9-i686-linux-gnu/bin/miranda {input.miR_Family_human_fasta} {params.output_directory}/temp.{wildcards.CHR}.fasta | gzip -c > {output.gencode_3utr_miRanda_output_for_a_single_chromosome}
        """

        
rule A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda____step03__get_edited_3UTR_sequences_from_edited_maf_blocks:
    input:
        edited_gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step07__get_edited_TargetScan_input_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/edited.gencode.3utr.{CHR}.maf.table.with.species.id.txt",
    output:
        A02_9__step02_flag=touch("result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step03__get_edited_3UTR_sequences_from_edited_maf_blocks/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished.{CHR}"),
        edited_gencode_3utr_fasta_for_a_single_chromosome="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step03__get_edited_3UTR_sequences_from_edited_maf_blocks/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/edited.gencode.3utr.{CHR}.fasta"
    shell:
        """
        bash ./scripts/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step03__get_edited_3UTR_sequences_from_edited_maf_blocks/convert_maf_to_human_fasta.sh {input.edited_gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename} > {output.edited_gencode_3utr_fasta_for_a_single_chromosome}
        """

rule A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda____step04__run_miRanda_per_edited_chromosome:
    input:
        edited_gencode_3utr_fasta_for_a_single_chromosome="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step03__get_edited_3UTR_sequences_from_edited_maf_blocks/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/edited.gencode.3utr.{CHR}.fasta",
        miR_Family_human_fasta="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step00__get_miRNA_Family_human_fasta/miR_Family_human.fasta"
    output:
        A02_9__step04_flag=touch("result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step04__run_miRanda_per_edited_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished.{CHR}"),
        edited_gencode_3utr_miRanda_output_for_a_single_chromosome="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step04__run_miRanda_per_edited_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/edited.gencode.3utr.{CHR}.miRanda.output.gz"
    params:
        output_directory="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step04__run_miRanda_per_edited_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/"
    shell:
        """
        mkdir -p {params.output_directory}
        rm -f {params.output_directory}/temp.{wildcards.CHR}.fasta
        ln -s -r {input.edited_gencode_3utr_fasta_for_a_single_chromosome} {params.output_directory}/temp.{wildcards.CHR}.fasta
        ./tools/miRanda-1.9-i686-linux-gnu/bin/miranda {input.miR_Family_human_fasta} {params.output_directory}/temp.{wildcards.CHR}.fasta | gzip -c > {output.edited_gencode_3utr_miRanda_output_for_a_single_chromosome}
        """


rule A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda____step05__concatenate_miRanda_results_across_all_chromosomes:
    input:
        miRanda_output_filenames_list=["result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step02__run_miRanda_per_chromosome/{GENCODE_VERSION}/gencode.3utr." + "chr" + str(chr) + ".miRanda.output.gz" for chr in list(range(1, 23)) + ["X", "Y"]]
    output:
        A02_9__step05_flag=touch("result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step05__concatenate_miRanda_results_across_all_chromosomes/{GENCODE_VERSION}/finished"),
        concatenated_headless_miRanda_output_gz_filename="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step05__concatenate_miRanda_results_across_all_chromosomes/{GENCODE_VERSION}/gencode.3utr.all.chromosomes.concatenated.headless.miRanda.output.gz"
    shell:
        """
        ## duplicated miRNA fasta names is OK; we only need the family info
        ## grep the "^>[^>]" lines, because they have 1-based [start, end] MBS on 3'-UTR (columns 8, 9)
        for file in {input.miRanda_output_filenames_list}
        do
            zcat $file | bash ./scripts/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step05__concatenate_miRanda_results_across_all_chromosomes/reformat.alignments.into.oneline.records.bash
        done | pigz -c > {output.concatenated_headless_miRanda_output_gz_filename}
        """


rule A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda____step06__concatenate_edited_miRanda_results_across_all_chromosomes:
    input:
        edited_miRanda_output_filenames_list=["result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step04__run_miRanda_per_edited_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/edited.gencode.3utr." + "chr" + str(chr) + ".miRanda.output.gz" for chr in list(range(1, 23)) + ["X"]]
    output:
        A02_9__step06_flag=touch("result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step06__concatenate_edited_miRanda_results_across_all_chromosomes/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished"),
        concatenated_headless_edited_miRanda_output_gz_filename="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step06__concatenate_edited_miRanda_results_across_all_chromosomes/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/edited.gencode.3utr.all.chromosomes.but.chrY.concatenated.headless.miRanda.output.gz"
    shell:
        """
        ## duplicated miRNA fasta names is OK; we only need the family info
        ## grep the "^>[^>]" lines, because they have 1-based [start, end] MBS on 3'-UTR (columns 8, 9)
        for file in {input.edited_miRanda_output_filenames_list}
        do
            zcat $file | bash ./scripts/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step06__concatenate_edited_miRanda_results_across_all_chromosomes/reformat.alignments.into.oneline.records.bash
        done | pigz -c > {output.concatenated_headless_edited_miRanda_output_gz_filename}
        """


rule A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda____step13__get_all_edited_3UTR_sequences_from_all_edited_maf_blocks:
    input:
        all_edited_gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename="result/A02_8__get_editing_effect_on_miRNA_binding_sites/step17__get_all_edited_TargetScan_input_per_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/all.edited.gencode.3utr.{CHR}.maf.table.with.species.id.txt",
    output:
        A02_9__step13_flag=touch("result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step13__get_all_edited_3UTR_sequences_from_edited_maf_blocks/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished.{CHR}"),
        all_edited_gencode_3utr_fasta_for_a_single_chromosome="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step13__get_all_edited_3UTR_sequences_from_edited_maf_blocks/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/all.edited.gencode.3utr.{CHR}.fasta"
    shell:
        """
        bash ./scripts/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step13__get_all_edited_3UTR_sequences_from_edited_maf_blocks/convert_maf_to_human_fasta.sh {input.all_edited_gencode_3utr_maf_table_with_species_id_for_a_single_chromosome_filename} > {output.all_edited_gencode_3utr_fasta_for_a_single_chromosome}
        """

rule A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda____step14__run_miRanda_per_all_edited_chromosome_per_miRNA:
    input:
        all_edited_gencode_3utr_fasta_for_a_single_chromosome="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step13__get_all_edited_3UTR_sequences_from_edited_maf_blocks/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/all.edited.gencode.3utr.{CHR}.fasta",
        miR_Family_human_fasta="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step00__get_miRNA_Family_human_fasta/miR_Family_human.fasta"
    output:
        A02_9__step14_flag=touch("result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step14__run_miRanda_per_all_edited_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/{MIRNA_NAME}/finished.{CHR}"),
        all_edited_gencode_3utr_miRanda_output_for_a_single_chromosome="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step14__run_miRanda_per_all_edited_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/{MIRNA_NAME}/all.edited.gencode.3utr.{CHR}.miRanda.output.gz"
    params:
        output_directory="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step14__run_miRanda_per_all_edited_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/{MIRNA_NAME}"
    shell:
        """
        mkdir -p {params.output_directory}
        rm -f {params.output_directory}/tools_for_{wildcards.CHR}
        ln -s -r ./tools {params.output_directory}/tools_for_{wildcards.CHR}
        rm -f {params.output_directory}/temp.miR_Family_human_for_{wildcards.CHR}.fasta
        grep -A 1 '>{wildcards.MIRNA_NAME}__' {input.miR_Family_human_fasta} > {params.output_directory}/temp.miR_Family_human_for_{wildcards.CHR}.fasta
        rm -f {params.output_directory}/temp.{wildcards.CHR}.fasta
        ln -s -r {input.all_edited_gencode_3utr_fasta_for_a_single_chromosome} {params.output_directory}/temp.{wildcards.CHR}.fasta
        cd {params.output_directory}
        ./tools_for_{wildcards.CHR}/miRanda-1.9-i686-linux-gnu/bin/miranda ./temp.miR_Family_human_for_{wildcards.CHR}.fasta ./temp.{wildcards.CHR}.fasta | gzip -c > ./all.edited.gencode.3utr.{wildcards.CHR}.miRanda.output.gz
        cd -
        """


miR_DataFrame = pandas.read_csv("./external/TargetScan/miR_Family_Info.txt", sep="\t")
list_of_miRNA_names = miR_DataFrame.loc[ (miR_DataFrame["Species ID"] == 9606) & (miR_DataFrame["Family Conservation?"] == 2), "MiRBase ID" ].tolist()

def write_input_list_to_file_and_return_this_filename(wildcards):
    all_edited_miRanda_output_filenames_list = ["result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step14__run_miRanda_per_all_edited_chromosome/" + wildcards.DATASET_COLLECTION_NAME + "/" + wildcards.DATASET_PHENOTYPE_COLLECTION_NAME + "/" + wildcards.SUBSET_NAME + "/" + wildcards.GENCODE_VERSION + "/" + miRNA_name + "/all.edited.gencode.3utr." + "chr" + str(chr) + ".miRanda.output.gz" for miRNA_name in list_of_miRNA_names for chr in list(range(1, 23)) + ["X", "Y"]]
    this_dir = "result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step16__concatenate_all_edited_miRanda_results_across_all_chromosomes/" + wildcards.DATASET_COLLECTION_NAME + "/" + wildcards.DATASET_PHENOTYPE_COLLECTION_NAME + "/" + wildcards.SUBSET_NAME + "/" + wildcards.GENCODE_VERSION
    os.makedirs(this_dir, exist_ok=True)
    this_file = this_dir + "/all.edited.miRanda.output.filenames"
    with open(this_file, "w") as f:
        [f.write(temp_filename + "\n") for temp_filename in all_edited_miRanda_output_filenames_list]
    return this_file
    
rule A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda____step16__concatenate_all_edited_miRanda_results_across_all_chromosomes:
    input:
        all_edited_miRanda_output_filenames_list=["result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step14__run_miRanda_per_all_edited_chromosome/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/" + miRNA_name + "/all.edited.gencode.3utr." + "chr" + str(chr) + ".miRanda.output.gz" for miRNA_name in list_of_miRNA_names for chr in list(range(1, 23)) + ["X", "Y"]]
    output:
        A02_9__step16_flag=touch("result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step16__concatenate_all_edited_miRanda_results_across_all_chromosomes/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/finished"),
        concatenated_headless_all_edited_miRanda_output_gz_filename="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step16__concatenate_all_edited_miRanda_results_across_all_chromosomes/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/all.edited.gencode.3utr.all.chromosomes.concatenated.headless.miRanda.output.gz"
    params:
        file_containing_all_edited_miRanda_output_filenames=write_input_list_to_file_and_return_this_filename,
        output_directory="result/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step16__concatenate_all_edited_miRanda_results_across_all_chromosomes/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{GENCODE_VERSION}/"
    shell:
        """
        ## duplicated miRNA fasta names is OK; we only need the family info
        ## grep the "^>[^>]" lines, because they have 1-based [start, end] MBS on 3'-UTR (columns 8, 9)
        mkdir -p {params.output_directory}
        rm -f {params.output_directory}/temp.out
        cat {params.file_containing_all_edited_miRanda_output_filenames} | while read file
        do
            echo `date` : processing $file
            zcat $file | bash ./scripts/A02_9__get_editing_effect_on_miRNA_binding_sites_for_miRanda/step16__concatenate_all_edited_miRanda_results_across_all_chromosomes/reformat.alignments.into.oneline.records.bash >> {params.output_directory}/temp.out
        done
        gzip -c {params.output_directory}/temp.out > {output.concatenated_headless_all_edited_miRanda_output_gz_filename}
        """
