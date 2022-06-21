for temp_key in ['threads_pigz']:
    if temp_key not in config.keys():
        config[temp_key] = 20

rule PS31_1__get_number_of_uniquely_mapped_bases:
    input:
        flag_PS15_1__step02__part05="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part05__recalibrate_base_quality",
    output:
        flag_PS31_1=touch("result/PS31_1__get_number_of_uniquely_mapped_bases/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished"),
        number_of_uniquely_mapped_bases="result/PS31_1__get_number_of_uniquely_mapped_bases/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/number.of.uniquely.mapped.bases"
    params:
        calling_result_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        working_result_directory="result/PS31_1__get_number_of_uniquely_mapped_bases/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/"
    shell:
        """
        ## same behavior for single-ended and paired-ended samples (for paired-ended samples, each NH:i:1 read pairs will have two and only two rows, and their total mapped bases is the very sum of the MD integers from these two rows)
        samtools view "{params.calling_result_directory}/alignment.sorted.withRG.dedup.reordered.recal.bam" | bash ./scripts/PS31_1__get_number_of_uniquely_mapped_bases/get_MDs_for_NH1.sh  | gzip -c > "{params.working_result_directory}/alignment.sorted.withRG.dedup.reordered.recal.bam.NH1.MDs.gz"
        zcat "{params.working_result_directory}/alignment.sorted.withRG.dedup.reordered.recal.bam.NH1.MDs.gz" | bash ./scripts/PS31_1__get_number_of_uniquely_mapped_bases/get_counts_from_MDs.sh > "{params.working_result_directory}/number.of.uniquely.mapped.bases"       
        """

rule PS32_1__get_basic_read_counts_for_mismatched_sites_only:
    input:
        flag_PS15_1__step02__part06="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part06__really_call_variants_in_rich_pileup_format",
        alignment_pileup_mismatched_sites_only_rich_filename="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/alignment.pileup.mismatched.sites.only.rich.gz"
    output:
        flag_PS32_1=touch("result/PS32_1__get_basic_read_counts_for_mismatched_sites_only/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{START_COORDINATE_OF_READ_3PRIME_END}/finished"),
        basic_read_counts_for_mismatched_sites_only_tsv_gz_filename="result/PS32_1__get_basic_read_counts_for_mismatched_sites_only/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{START_COORDINATE_OF_READ_3PRIME_END}/basic.read.counts.for.mismatched.sites.only.tsv.gz"
    threads:
        config['threads_pigz']
    shell:
        """
        pigz -p {threads} -d -c {input.alignment_pileup_mismatched_sites_only_rich_filename} | cut -f 1-6,8 | python ./scripts/PS32_1__get_basic_read_counts_for_mismatched_sites_only/filter.1.2.3.get.read.counts.py {wildcards.START_COORDINATE_OF_READ_3PRIME_END} | pigz -c > {output.basic_read_counts_for_mismatched_sites_only_tsv_gz_filename}
        """

rule PS33_1__get_base_and_read_name_pairs_for_mismatched_sites_only:
    input:
        flag_PS15_1__step02__part06="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part06__really_call_variants_in_rich_pileup_format",
        alignment_pileup_mismatched_sites_only_rich_filename="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/alignment.pileup.mismatched.sites.only.rich.gz"
    output:
        flag_PS33_1=touch("result/PS33_1__get_base_and_read_name_pairs_for_mismatched_sites_only/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished"),
        base_and_read_name_pairs_for_mismatched_sites_only_tsv_gz_filename="result/PS33_1__get_base_and_read_name_pairs_for_mismatched_sites_only/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/base.and.read.name.pairs.for.mismatched.sites.only.tsv.gz"
    threads:
        config['threads_pigz']
    shell:
        """
        pigz -p {threads} -d -c {input.alignment_pileup_mismatched_sites_only_rich_filename} | cut -f 1-4,5,7 | python ./scripts/PS33_1__get_base_and_read_name_pairs_for_mismatched_sites_only/filter.4.get.variant.supporting.reads.py | pigz -c > {output.base_and_read_name_pairs_for_mismatched_sites_only_tsv_gz_filename}
        """
