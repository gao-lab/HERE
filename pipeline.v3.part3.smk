import pandas

for temp_key in ['threads_reduce_and_pigz_compress_tables', 'threads_concatenating_vcfs', 'threads_bcftools_isec', 'threads_filter_for_variants_with_enough_read_support', 'threads_filter_for_variants_with_enough_sample_support', 'threads_filter_for_A_to_G_only_variants', 'threads_mark_unsequenced_editing_sites', 'threads_annotate_embryonic_genes']:
    if temp_key not in config.keys():
        config[temp_key] = 1


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


## must take unique at the depth filename list step to avoid duplicated records
def __S41_1_step03_write_alignment_merged_depth_on_merged_variant_only_bed_txt_gz_filename_collection_to_file_and_return_filename(wildcards):
    temp_result_directory = "result/S41_1__check_variant_converage_of_merged_bam/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + wildcards.DATASET_COLLECTION_NAME_FOR_READ_COVERAGE_CHECK + "/"
    temp_dataset_and_sample_and_depth_filename_list_filename = temp_result_directory + "/temp_dataset_and_sample_and_depth_filename_list.txt"
    os.makedirs(temp_result_directory, exist_ok=True)
    temp_dataset_and_sample_and_depth_filenames_collection = set([temp_row.DATASET_NAME + "," + temp_row.SAMPLE_NAME + "," + "result/S41_1__check_variant_converage_of_merged_bam/" + wildcards.DATASET_COLLECTION_NAME + "/__merged__/" + wildcards.FASTP_PARAMETERS + "/" + wildcards.CONTIGS_FASTA_FILENAME + "/" + wildcards.GENCODE_VERSION + "/" + wildcards.INDEXER + "/" + wildcards.ALIGNER + "/" + wildcards.ALIGNER_PARAMETERS + "/" + wildcards.VARIANT_CALLER + "/" + wildcards.VARIANT_CALLER_PARAMETERS + "/" + wildcards.dbSNP_VERSION + "/" + wildcards.dbSNP_SUBSET + "/" + wildcards.COMPLEX_FILTER + "/" + wildcards.COMPLEX_FILTER_PARAMETERS + "/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + "/" + temp_row.SAMPLE_NAME + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/alignment.merged.depth.on.merged.variant.only.bed.txt.gz" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME_FOR_READ_COVERAGE_CHECK)]) ## take unique datasets
    with open(temp_dataset_and_sample_and_depth_filename_list_filename, "w") as temp_self_filehandle:
        print("writing...")
        temp_self_filehandle.writelines([temp_dataset_and_sample_and_depth_filename + "\n" for temp_dataset_and_sample_and_depth_filename in temp_dataset_and_sample_and_depth_filenames_collection])
    return temp_dataset_and_sample_and_depth_filename_list_filename



rule S42_1__annotate_embryonic_genes:
    input:
        flag_S06_1__step03="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/finished.step03__get_expression_matrix_by_ballgown",
        combined_gexpr_FPKM_matrix_filename="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/combined.gexpr.FPKM.matrix.txt",
        flag_S21_1="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        phenotype_output_at_gsm_level_dt_filename="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/phenotype.output.at.gsm.level.dt.txt",
        reference_GTF_filename="external/reference.gene.annotation/GENCODE.annotation/{GENCODE_VERSION}/gencode.annotation.gtf"
    output:
        flag_S42_1=touch("result/S42_1__annotate_embryonic_genes/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished"),
        combined_gexpr_FPKM_pc_only_melt_with_phenotype_normal_sample_only_dt_txt_gz_filename="result/S42_1__annotate_embryonic_genes/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/{DATASET_PHENOTYPE_COLLECTION_NAME}/combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.dt.txt.gz",
        combined_gexpr_FPKM_pc_only_melt_with_phenotype_normal_sample_only_median_annotated_dt_txt_gz_filename="result/S42_1__annotate_embryonic_genes/{DATASET_COLLECTION_NAME}/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/__sample_dependent__/{ALIGNER}/{ALIGNER_PARAMETERS}/{EXPRESSION_CALLER}/{EXPRESSION_CALLER_PARAMETERS}/{DATASET_PHENOTYPE_COLLECTION_NAME}/combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt.txt.gz"
    threads:
        config['threads_annotate_embryonic_genes']
    script:
        "./scripts/S42_1__annotate_embryonic_genes/run.R"





rule S51_1__combine_multiple_population_isec_results:
    input:
        isec_UWashington_EVS_vcf_gz_filename="result/S19_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/bcftools.isec.with.outer.vcf/basic/10000000/UWashington.EVS/collapse_all_and_keep_self_vcf/combined.merged.variant.only.bcftools.isec.outer.vcf.result.self.vcf.gz",
        isec_NCBI_ALFA_vcf_gz_filename="result/S19_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/bcftools.isec.with.outer.vcf/basic/10000000/NCBI.ALFA.2020.03.04/collapse_all_and_keep_self_vcf/combined.merged.variant.only.bcftools.isec.outer.vcf.result.self.vcf.gz",
        isec_1000Genomes_vcf_gz_filename="result/S19_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/bcftools.isec.with.outer.vcf/basic_without_Y/10000000/1000Genomes.phased.genotype/collapse_all_and_keep_self_vcf/combined.merged.variant.only.bcftools.isec.outer.vcf.result.self.vcf.gz",
        isec_gnomAD_v211_exomes_vcf_gz_filename="result/S19_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/bcftools.isec.with.outer.vcf/basic/10000000/gnomAD_v2.1.1_exomes/collapse_all_and_keep_self_vcf/combined.merged.variant.only.bcftools.isec.outer.vcf.result.self.vcf.gz",
        isec_gnomAD_v211_genomes_vcf_gz_filename="result/S19_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/bcftools.isec.with.outer.vcf/basic_without_Y/10000000/gnomAD_v2.1.1_genomes/collapse_all_and_keep_self_vcf/combined.merged.variant.only.bcftools.isec.outer.vcf.result.self.vcf.gz",
        isec_gnomAD_v30_genomes_vcf_gz_filename="result/S19_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/bcftools.isec.with.outer.vcf/basic/10000000/gnomAD_v3.0_genomes/collapse_all_and_keep_self_vcf/combined.merged.variant.only.bcftools.isec.outer.vcf.result.self.vcf.gz",
    output:
        flag_S51_1=touch("result/S51_1__combine_multiple_population_isec_results/{DATASET_COLLECTION_NAME}/finished"),
        isec_all_populations_concatenated_vcf_gz_filename="result/S51_1__combine_multiple_population_isec_results/{DATASET_COLLECTION_NAME}/isec.all.populations.concatenated.vcf.gz"
    threads:
        config['threads_concatenating_vcfs']
    shell:
        """
        bcftools concat --threads {threads} --allow-overlaps -o {output.isec_all_populations_concatenated_vcf_gz_filename} -Oz {input.isec_UWashington_EVS_vcf_gz_filename} {input.isec_NCBI_ALFA_vcf_gz_filename} {input.isec_1000Genomes_vcf_gz_filename} {input.isec_gnomAD_v211_exomes_vcf_gz_filename} {input.isec_gnomAD_v211_genomes_vcf_gz_filename} {input.isec_gnomAD_v30_genomes_vcf_gz_filename}
        tabix {output.isec_all_populations_concatenated_vcf_gz_filename}
        """

rule S51_2__filter_against_population_variants:
    input:
        merged_variant_only_vcf_gz_filename="result/S16_1__concatenate_RNA_editing_site_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/merged.variant.only.vcf.gz",
        isec_all_populations_concatenated_vcf_gz_filename="result/S51_1__combine_multiple_population_isec_results/{DATASET_COLLECTION_NAME}/isec.all.populations.concatenated.vcf.gz"
    output:
        flag_S51_2=touch("result/S51_2__filter_against_population_variants/{DATASET_COLLECTION_NAME}/finished"),
        merged_variant_only_disjoint_with_population_variants_vcf_gz_filename=touch("result/S51_2__filter_against_population_variants/{DATASET_COLLECTION_NAME}/merged.variant.only.disjoint.with.population.variants.vcf.gz")
    params:
        temp_result_directory="result/S51_2__filter_against_population_variants/{DATASET_COLLECTION_NAME}/temp_result"
    threads:
        config['threads_bcftools_isec']
    shell:
        """
        rm -fr {params.temp_result_directory}/*
        mkdir -p {params.temp_result_directory}
        bcftools isec --threads {threads} --collapse all --output-type z --prefix {params.temp_result_directory}/ {input.merged_variant_only_vcf_gz_filename} {input.isec_all_populations_concatenated_vcf_gz_filename}
        ln -s -r {params.temp_result_directory}/0000.vcf.gz {output.merged_variant_only_disjoint_with_population_variants_vcf_gz_filename}
        tabix {output.merged_variant_only_disjoint_with_population_variants_vcf_gz_filename}
        """

rule S51_3__filter_for_variants_with_enough_read_support:
    input:
        flag_S51_2="result/S51_2__filter_against_population_variants/{DATASET_COLLECTION_NAME}/finished",
        merged_variant_only_disjoint_with_population_variants_vcf_gz_filename="result/S51_2__filter_against_population_variants/{DATASET_COLLECTION_NAME}/merged.variant.only.disjoint.with.population.variants.vcf.gz",
        merged_long_table_dt_txt_gz_filename="result/S16_3__get_RNA_editing_site_long_table_from_a_dataset_collection/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/merged.long.table.dt.txt.gz",
        combined_merged_variant_only_snpEff_event_summary_dt_gz_filename="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/snpEff/basic/10000000/combined.merged.variant.only.snpEff.event.summary.dt.txt.gz"
    output:
        flag_S51_3=touch("result/S51_3__filter_for_variants_with_enough_read_support/{DATASET_COLLECTION_NAME}/finished"),
        merged_long_disjoint_with_population_dt_txt_gz_filename="result/S51_3__filter_for_variants_with_enough_read_support/{DATASET_COLLECTION_NAME}/merged.long.disjoint.with.population.dt.txt.gz",
        merged_long_disjoint_with_population_without_potential_polymorphism_dt_txt_gz_filename="result/S51_3__filter_for_variants_with_enough_read_support/{DATASET_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.dt.txt.gz",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_dt_txt_gz_filename="result/S51_3__filter_for_variants_with_enough_read_support/{DATASET_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt.txt.gz"
    threads:
        config['threads_filter_for_variants_with_enough_read_support']
    script:
        "./scripts/S51_3__filter_for_variants_with_enough_read_support/run.R"




rule S51_4__filter_for_variants_with_enough_sample_support:
    input:
        flag_S51_3="result/S51_3__filter_for_variants_with_enough_read_support/{DATASET_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_dt_txt_gz_filename="result/S51_3__filter_for_variants_with_enough_read_support/{DATASET_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.dt.txt.gz",
        flag_S21_1="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        phenotype_output_at_gsm_level_dt_filename="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/phenotype.output.at.gsm.level.dt.txt"
    output:
        flag_S51_4=touch("result/S51_4__filter_for_variants_with_enough_sample_support/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished"),
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_dt_txt_gz_filename="result/S51_4__filter_for_variants_with_enough_sample_support/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.dt.txt.gz",
        phenotype_of_samples_without_edits_before_sample_support_filter_dt_filename="result/S51_4__filter_for_variants_with_enough_sample_support/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/phenotype.of.samples.without.edits.before.sample.support.filter.dt.txt",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_dt_txt_gz_filename="result/S51_4__filter_for_variants_with_enough_sample_support/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt.txt.gz"
    threads:
        config['threads_filter_for_variants_with_enough_sample_support']
    script:
        "./scripts/S51_4__filter_for_variants_with_enough_sample_support/run.R"

rule S51_5__filter_for_A_to_G_sites:
    input:
        flag_S51_4="result/S51_4__filter_for_variants_with_enough_sample_support/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_dt_txt_gz_filename="result/S51_4__filter_for_variants_with_enough_sample_support/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.dt.txt.gz",
        combined_merged_variant_only_snpEff_event_summary_dt_gz_filename="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/snpEff/basic/10000000/combined.merged.variant.only.snpEff.event.summary.dt.txt.gz"
    output:
        flag_S51_5=touch("result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished"),
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_with_event_summary_dt_txt_gz_filename="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_bed_filename="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.bed"
    threads:
        config['threads_filter_for_A_to_G_only_variants']
    script:
        "./scripts/S51_5__filter_for_A_to_G_sites/run.R"




rule S51_6__get_snpEff_annotation_subset_of_filtered_result:
    input:
        flag_S51_5="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_with_event_summary_dt_txt_gz_filename="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz",
        ## borrow the whole SnpEff annotation result
        flag_S18_1_step02_fixed="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/snpEff/basic/10000000/finished.step02__combine_merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_filename____patch01__get_full_annotation",
        combined_merged_vcf_reformatted_with_snpEff_ANN_split_annotation_dt_gz_filename="result/S18_1__combine_annotations/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/complex_filter_1/none/snpEff/basic/10000000/combined.merged.variant.only.snpEff.ANN.single.match.dt.txt.gz"
    output:
        flag_S51_6=touch("result/S51_6__get_snpEff_annotation_subset_of_filtered_result/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished"),
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_with_snpEff_annotation_dt_txt_gz_filename="result/S51_6__get_snpEff_annotation_subset_of_filtered_result/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.with.snpEff.annotation.dt.txt.gz",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_with_event_summary_variant_only_snpEff_annotation_dt_txt_gz_filename="result/S51_6__get_snpEff_annotation_subset_of_filtered_result/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.variant.only.snpEff.annotation.dt.txt.gz",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_variant_only_snpEff_annotation_dt_txt_gz_filename="result/S51_6__get_snpEff_annotation_subset_of_filtered_result/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt.txt.gz"
    script:
        "./scripts/S51_6__get_snpEff_annotation_subset_of_filtered_result/run.R"




rule S52_1__check_variant_converage_of_merged_bam:
    input:
        flag_S51_5="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_bed_filename="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.bed",
        alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/{INDEXER_PARAMETERS}/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
    output:
        flag_S52_1=touch("result/S52_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished"),
        alignment_merged_depth_on_merged_variant_only_bed_txt_gz_filename="result/S52_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.merged.variant.only.bed.txt.gz"
    threads:
        1
    shell:
        """
        samtools depth -a -d 0 -b {input.merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_bed_filename} {input.alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename} | gzip > {output.alignment_merged_depth_on_merged_variant_only_bed_txt_gz_filename}
        """


rule B52_1__check_variant_converage_of_merged_bam:
    input:
        lambda wildcards: ["result/S52_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/finished" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_B52_1=touch("result/B52_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished")




## must take unique at the depth filename list step to avoid duplicated records
def __S52_2_write_alignment_merged_depth_on_merged_variant_only_bed_txt_gz_filename_collection_to_file_and_return_filename(wildcards):
    temp_result_directory = "result/S52_2__concatenate_all_variant_coverages_of_merged_bam/" + wildcards.DATASET_COLLECTION_NAME + "/"  + wildcards.DATASET_PHENOTYPE_COLLECTION_NAME  + "/"
    temp_dataset_and_sample_and_depth_filename_list_filename = temp_result_directory + "/temp_dataset_and_sample_and_depth_filename_list.txt"
    os.makedirs(temp_result_directory, exist_ok=True)
    temp_dataset_and_sample_and_depth_filenames_collection = set([temp_row.DATASET_NAME + "," + temp_row.SAMPLE_NAME + "," + "result/S52_1__check_variant_converage_of_merged_bam/" + wildcards.DATASET_COLLECTION_NAME + "/"  + wildcards.DATASET_PHENOTYPE_COLLECTION_NAME  + "/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + "/" + temp_row.SAMPLE_NAME + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/alignment.merged.depth.on.merged.variant.only.bed.txt.gz" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]) ## take unique datasets
    with open(temp_dataset_and_sample_and_depth_filename_list_filename, "w") as temp_self_filehandle:
        print("writing...")
        temp_self_filehandle.writelines([temp_dataset_and_sample_and_depth_filename + "\n" for temp_dataset_and_sample_and_depth_filename in temp_dataset_and_sample_and_depth_filenames_collection])
    return temp_dataset_and_sample_and_depth_filename_list_filename


rule S52_2__concatenate_all_variant_coverages_of_merged_bam:
    input:
        flag_B52_1="result/B52_1__check_variant_converage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished"
    output:
        flag_S52_2=touch("result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished"),
        combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_gz_filename="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.txt.gz"
    params:
        alignment_merged_depth_on_merged_variant_only_bed_txt_gz_filename_collection_filename=lambda wildcards: __S52_2_write_alignment_merged_depth_on_merged_variant_only_bed_txt_gz_filename_collection_to_file_and_return_filename(wildcards),
        temp_combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_filename="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.txt"
    threads:
        config['threads_reduce_and_pigz_compress_tables']
    shell:
        """
        rm -f {params.temp_combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_filename}
        touch {params.temp_combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_filename}
        cat {params.alignment_merged_depth_on_merged_variant_only_bed_txt_gz_filename_collection_filename} | while read line; do
            ## dataset=`echo $line|cut -f 1 -d ,`
            sample=`echo $line|cut -f 2 -d ,`
            depth_filename=`echo $line|cut -f 3 -d ,`
            zcat $depth_filename | awk -v sample=$sample '{{print sample"\t"$0}}' >> {params.temp_combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_filename}
        done
        pigz --processes {threads} --stdout {params.temp_combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_filename} > {output.combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_gz_filename}
        """



rule S52_2__concatenate_all_variant_coverages_of_merged_bam____patch01__extract_zero_depth_records:
    input:
        flag_S52_2="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_gz_filename="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.txt.gz"
    output:
        flag_S52_2_patch01=touch("result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished.patch01__extract_zero_depth_records"),
        combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_zero_depth_records_only_gz_filename="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.zero.depth.records.only.txt.gz"
    threads:
        config['threads_reduce_and_pigz_compress_tables']
    shell:
        """
        zcat {input.combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_gz_filename} | awk '{{if($4 == 0) print}}' | pigz --processes {threads} --stdout > {output.combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_zero_depth_records_only_gz_filename}
        """


rule S52_2__concatenate_all_variant_coverages_of_merged_bam____patch02__extract_nonzero_depth_records:
    input:
        flag_S52_2="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_gz_filename="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.txt.gz"
    output:
        flag_S52_2_patch02=touch("result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished.patch02__extract_nonzero_depth_records"),
        combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_nonzero_depth_records_only_gz_filename="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.nonzero.depth.records.only.txt.gz"
    threads:
        config['threads_reduce_and_pigz_compress_tables']
    shell:
        """
        zcat {input.combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_gz_filename} | awk '{{if($4 > 0) print}}' | pigz --processes {threads} --stdout > {output.combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_txt_nonzero_depth_records_only_gz_filename}
        """


rule S52_3__mark_unsequenced_editing_sites:
    input:
        flag_S51_5="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz",
        flag_S52_2_patch01="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished.patch01__extract_zero_depth_records",
        combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_zero_depth_records_only_txt_gz_filename="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.zero.depth.records.only.txt.gz",
        flag_S52_2_patch02="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished.patch02__extract_nonzero_depth_records",
        combined_all_single_bam_coverages_of_the_bed_generated_from_vcf_nonzero_depth_records_only_txt_gz_filename="result/S52_2__concatenate_all_variant_coverages_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/combined.all.single.bam.coverages.of.the.bed.generated.from.vcf.nonzero.depth.records.only.txt.gz",
        flag_S21_1="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        phenotype_output_at_gsm_level_dt_filename="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/phenotype.output.at.gsm.level.dt.txt"
    output:
        flag_S52_3=touch("result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished"),
        test_merge_result_dt_txt_filename="result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/test.merge.result.dt.txt",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_merged_with_coverage_dt_txt_gz_filename="result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt.txt.gz"
    threads:
        config['threads_mark_unsequenced_editing_sites']
    script:
        "./scripts/S52_3__mark_unsequenced_editing_sites/run.R"


rule S53_1__generate_expression_summary_with_phenotype:
    input:
        flag_S06_1__step03="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/finished.step03__get_expression_matrix_by_ballgown",
        combined_gexpr_FPKM_matrix_filename="result/BS06_1__get_expression_level/{DATASET_COLLECTION_NAME}/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/combined.gexpr.FPKM.matrix.txt",
        flag_S21_1="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        phenotype_output_at_gsm_level_dt_filename="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/phenotype.output.at.gsm.level.dt.txt"
    output:
        flag_S53_1=touch("result/S53_1__generate_expression_summary_with_phenotype/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished"),
        aa="result/S53_1__generate_expression_summary_with_phenotype/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/??"
    script:
        "./scripts/S53_1__generate_expression_summary_with_phenotype/run.R"
