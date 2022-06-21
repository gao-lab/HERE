for temp_key in ['threads_motif_calling', 'threads_auxiliary_processing']:
    if temp_key not in config.keys():
        config[temp_key] = 1


rule A01_1__compute_total_sample_count:
    input:
        flag_S21_1="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        phenotype_output_at_gsm_level_dt_filename="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/phenotype.output.at.gsm.level.dt.txt",
    output:
        flag=touch("report.ver2/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/total.sample.count.for.normal.stages.finished"),
        total_sample_count_for_normal_stages_dt_csv_filename="report.ver2/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/total.sample.count.for.normal.stages.dt.csv",
        total_sample_count_for_GSE133854_all_dt_csv_filename="report.ver2/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/total.sample.count.for.GSE133854.all.dt.csv"
    threads:
        1
    script:
        "./scripts/A01_1__compute_total_sample_count/run.R"


rule A01_5__plot_motif:
    input:
        flag_S52_3="result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_bed_filename="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.bed",
        reference_genome_fasta_fai_filename="external/contigs/hg38.fa.fai",
        reference_genome_fasta_filename="external/contigs/hg38.fa",
        reference_transcriptome_fasta_gz_filename="external/contigs/gencode.v32.transcripts.fa.gz"
    output:
        flag_A01_5=touch("result/A01_5__plot_motif/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished")
    params:
        temp_result_directory="result/A01_5__plot_motif/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/temp_result_directory"
    threads:
        config["threads_motif_calling"]
    shell:
        """
        rm -fr {params.temp_result_directory}/*
        mkdir -p {params.temp_result_directory}
        cat {input.merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_bed_filename} | awk 'BEGIN{{OFS="\t"}} {{if ($4 ~ /A_G/) {{$6="+"}} else {{$6="-"}}; print }}' | sort -k1.4,2 -V  > {params.temp_result_directory}/temp.stranded.sorted.bed
        bedtools slop -i {params.temp_result_directory}/temp.stranded.sorted.bed -g {input.reference_genome_fasta_fai_filename} -b 3 > {params.temp_result_directory}/temp.b3.extended.bed
        bedtools getfasta -s -fi {input.reference_genome_fasta_filename} -bed {params.temp_result_directory}/temp.b3.extended.bed -fo {params.temp_result_directory}/temp.b3.extended.fasta
        grep -v '>' {params.temp_result_directory}/temp.b3.extended.fasta | tr [a-z] [A-Z] > {params.temp_result_directory}/temp.b3.extended.sequence.lines

        seqkit locate --only-positive-strand --use-regexp --pattern "...A..." --threads {threads} {input.reference_transcriptome_fasta_gz_filename} | tail -n +2 | cut -f 7 > {params.temp_result_directory}/background.7mer.lines

        cp -r tools/tsl {params.temp_result_directory}/
        cd {params.temp_result_directory}/tsl/cgi-bin/
        ./tsl -P ../../temp.b3.extended.sequence.lines -N ../../background.7mer.lines -K N -O ../../temp.b3.extended.sequence.lines.tsl.plot.png -x -y -C nucleo_weblogo -U cm -W 8 -H 4 -R 300
        cd -
        """
