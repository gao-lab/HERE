rule PS80_1__index_contig_with_bwa:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}"
    output:
        flag_PS80_1=touch("result/PS80_1__index_contig_with_bwa/{CONTIGS_FASTA_FILENAME}/{INDEXER}/{INDEXER_PARAMETERS}/finished")
    params:
        contigs_index_directory="result/PS80_1__index_contig_with_bwa/{CONTIGS_FASTA_FILENAME}/{INDEXER}/{INDEXER_PARAMETERS}/",
        contigs_index_prefix="result/PS80_1__index_contig_with_bwa/{CONTIGS_FASTA_FILENAME}/{INDEXER}/{INDEXER_PARAMETERS}/prefix"
    threads:
        1
    shell:
        """
        rm -fr "{params.contigs_index_directory}"/*
        mkdir -p "{params.contigs_index_directory}"
        if [[ "{wildcards.INDEXER_PARAMETERS}" == 'bwa-index-default' ]]; then
            ./tools/bwa-0.6.2/bwa index -p {params.contigs_index_prefix} {input.contigs_fasta_filename}
        else
            echo "Indexer parameters {wildcards.INDEXER_PARAMETERS} not supported." && exit 1;
        fi
        """


rule PS80_2__generate_snpdb:
    input:
        snpdb="external/snpdb/{SNPDB}"
    output:
        flag_PS80_2=touch("result/PS80_2__generate_snpdb/{SNPDB}/finished"),
        snp_temp_txt="result/PS80_2__generate_snpdb/{SNPDB}/snp.temp.txt"
    params:
        output_directory="result/PS80_2__generate_snpdb/{SNPDB}"
    threads:
        1
    shell:
        """
        perl scripts/PS80_2__generate_snpdb/build_snpdb.pl --snpdb {input.snpdb} --output {params.output_directory}
        """



rule PS81_1__Qiu2016_call_variants:
    input:
        flag_PS15_1__step02__part05="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part05__recalibrate_base_quality",
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}"
    output:
        flag_PS81_1=touch("result/PS81_1__Qiu2016_call_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished")
    params:
        bam_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        call_directory="result/PS81_1__Qiu2016_call_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/"
    threads:
        2
    shell:
        """
        ./external/Qiu2016.rep/RNA_editing_pipeline/bin/MismatchStat -i "{params.bam_directory}"/alignment.sorted.withRG.dedup.reordered.recal.bam -x 1000000 -u -o "{params.call_directory}"/rna_mismatch_stat.txt
        ./external/Qiu2016.rep/RNA_editing_pipeline/bin/MutDet -i "{params.bam_directory}"/alignment.sorted.withRG.dedup.reordered.recal.bam -r {input.contigs_fasta_filename}  -q 50 -v "{params.call_directory}"/rna_mismatch_stat.txt -u -o "{params.call_directory}"/variation.sites.txt.gz
        """




rule PS81_2__Qiu2016_filter_variants____part1__binom_and_readend_and_strand_and_mism_and_poly_and_rep_and_snp:
    input:
        flag_PS15_1__step02__part05="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished.step02__call_variants____part05__recalibrate_base_quality",
        flag_PS81_1="result/PS81_1__Qiu2016_call_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/finished",
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}",
        hg19_simple_repeat_bed="external/UCSC.Table.Browser/hg19_simpleRepeat.reg.bed",
        hg38_simple_repeat_bed="external/UCSC.Table.Browser/hg38_simpleRepeat.reg.bed",
        flag_PS80_2="result/PS80_2__generate_snpdb/{SNPDB}/finished",
        snp_temp_txt="result/PS80_2__generate_snpdb/{SNPDB}/snp.temp.txt"
    output:
        flag_PS81_2__part1=touch("result/PS81_2__Qiu2016_filter_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{SNPDB}/finished.part1__binom_and_readend_and_strand_and_mism_and_poly_and_rep_and_snp")
    params:
        bam_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        call_directory="result/PS81_1__Qiu2016_call_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        filter_directory="result/PS81_2__Qiu2016_filter_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{SNPDB}/",
        perl_scripts_directory="external/Qiu2016.rep/RNA_editing_pipeline/bin/"
    threads:
        1
    shell:
        """
        hgassembly_simple_repeat_bed=""
        if [[ {wildcards.CONTIGS_FASTA_FILENAME} == "hg19.fa" ]] ; then
            hgassembly_simple_repeat_bed={input.hg19_simple_repeat_bed}
        elif [[ {wildcards.CONTIGS_FASTA_FILENAME} == "hg38.fa" ]] ; then
            hgassembly_simple_repeat_bed={input.hg38_simple_repeat_bed}
        else 
            echo "Unsupported CONTIGS_FASTA_FILENAME: {wildcards.CONTIGS_FASTA_FILENAME}" && exit 1
        fi
        

        #Binomial,variation reads,editing frequency,strandbias,reads end,simple repeat and homopolymer filter;
        perl {params.perl_scripts_directory}/binom.reads.fre.end.strand.mism.poly.rep.filter.pl --in {params.call_directory}/variation.sites.txt.gz --out {params.filter_directory}/binom.reads.fre.strand.end.simplerepeat.homopolymer.txt --simpleRepeat $hgassembly_simple_repeat_bed --reference {input.contigs_fasta_filename} --endfre 0.9 --endp 0.05 --strandfre 0.9 --strandp 0.005 --homopolymer 5
        
        #dbSNP filter;
        perl {params.perl_scripts_directory}/snp.filter.pl {params.filter_directory}/binom.reads.fre.strand.end.simplerepeat.homopolymer.txt {input.snp_temp_txt} {params.filter_directory}/binom.reads.fre.strand.end.simplerepeat.homopolymer.dbsnp.txt

        ln -s -r {params.filter_directory}/binom.reads.fre.strand.end.simplerepeat.homopolymer.dbsnp.txt {params.filter_directory}/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.txt
        """


rule PS81_2__Qiu2016_filter_variants____part2__bwa_realignment_and_filtering:
    input:
        flag_PS81_2__part1="result/PS81_2__Qiu2016_filter_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{SNPDB}/finished.part1__binom_and_readend_and_strand_and_mism_and_poly_and_rep_and_snp",
        flag_PS80_1_hg19="result/PS80_1__index_contig_with_bwa/hg19.fa.and.Homo_sapiens.GRCh37.75.cdna.all.fa/bwa/bwa-index-default/finished",
        flag_PS80_1_hg38="result/PS80_1__index_contig_with_bwa/hg38.fa.and.gencode.v32.transcripts.fa/bwa/bwa-index-default/finished",
    output:
        flag_PS81_2__part2=touch("result/PS81_2__Qiu2016_filter_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{SNPDB}/finished.part2__bwa_realignment_and_filtering")
    params:
        bam_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        call_directory="result/PS81_1__Qiu2016_call_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        filter_directory="result/PS81_2__Qiu2016_filter_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{SNPDB}",
        perl_scripts_directory="external/Qiu2016.rep/RNA_editing_pipeline/bin/",
        hg19_mRNA_contigs_index_prefix="result/PS80_1__index_contig_with_bwa/hg19.fa.and.Homo_sapiens.GRCh37.75.cdna.all.fa/bwa/bwa-index-default/prefix",
        hg38_mRNA_contigs_index_prefix="result/PS80_1__index_contig_with_bwa/hg38.fa.and.gencode.v32.transcripts.fa/bwa/bwa-index-default/prefix",
        bwa_sampe_option_a=170*1.3,
        annovar_directory="tools/annovar/"
    threads:
        1
    shell:
        """
        Bin="{params.perl_scripts_directory}"
        rnabam="{params.bam_directory}"/alignment.sorted.withRG.dedup.reordered.recal.bam
        output={params.filter_directory}

        hgassembly_mRNA=""
        annovar_humandb_refGene=""
        if [[ {wildcards.CONTIGS_FASTA_FILENAME} == "hg19.fa" ]] ; then
            hgassembly_mRNA={params.hg19_mRNA_contigs_index_prefix}
            annovar_humandb_refGene={params.annovar_directory}/humandb_hg19/hg19_refGene.txt
        elif [[ {wildcards.CONTIGS_FASTA_FILENAME} == "hg38.fa" ]] ; then
            hgassembly_mRNA={params.hg38_mRNA_contigs_index_prefix}
            annovar_humandb_refGene={params.annovar_directory}/humandb_hg38/hg38_refGene.txt 
        else 
            echo "Unsupported CONTIGS_FASTA_FILENAME: {wildcards.CONTIGS_FASTA_FILENAME}" && exit 1
        fi

        type_end=`echo {wildcards.TYPE} | cut -f 1 -d '-'`
        
        ## Extract reads in which candidate RNA editing sites from the rna bam...
        mkdir -p $output/bwa


        if [[ $type_end == 'paired' ]]; then
            perl $Bin/extract_reads.pl $rnabam $output/bwa/mutation.read1.fq $output/bwa/mutation.read2.fq $output/bwa/mutation.read.fq $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.txt $output/readsvnum.txt
        elif [[ $type_end == 'single' ]]; then
            perl $Bin/extract_reads_SE-1.pl $rnabam $output/bwa/mutation.read.fq $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.txt $output/readsvnum.txt
        else
            echo "Unsupported type: {wildcards.TYPE}" && exit 1       
        fi

        ## Realign the reads to the combined reference...
        if [[ $type_end == 'paired' ]]; then
            tools/bwa-0.6.2/bwa aln -o 1 -l 31 -k 2 -t 4 -L -i 15 $hgassembly_mRNA $output/bwa/mutation.read1.fq > $output/bwa/mutation.read1.sai
            tools/bwa-0.6.2/bwa aln -o 1 -l 31 -k 2 -t 4 -L -i 15 $hgassembly_mRNA $output/bwa/mutation.read2.fq > $output/bwa/mutation.read2.sai
            tools/bwa-0.6.2/bwa sampe -a {params.bwa_sampe_option_a} $hgassembly_mRNA $output/bwa/mutation.read1.sai $output/bwa/mutation.read2.sai $output/bwa/mutation.read1.fq $output/bwa/mutation.read2.fq | samtools view -b -S - -t $hgassembly_mRNA.fai > $output/bwa/mutation.readP.bam
            tools/bwa-0.6.2/bwa aln -o 1 -e 50 -m 100000 -l 32 -k 2 -t 4 -L -i 15 $hgassembly_mRNA $output/bwa/mutation.read.fq > $output/bwa/mutation.read.sai
            tools/bwa-0.6.2/bwa samse $hgassembly_mRNA $output/bwa/mutation.read.sai $output/bwa/mutation.read.fq  | samtools view -b -S -t $hgassembly_mRNA.fai - > $output/bwa/mutation.readS.bam
        elif [[ $type_end == 'single' ]]; then
            tools/bwa-0.6.2/bwa aln -o 1 -e 50 -m 100000 -l 32 -k 2 -t 4 -L -i 15 $hgassembly_mRNA $output/bwa/mutation.read.fq > $output/bwa/mutation.read.sai
            tools/bwa-0.6.2/bwa samse $hgassembly_mRNA $output/bwa/mutation.read.sai $output/bwa/mutation.read.fq  | samtools view -b -S -t $hgassembly_mRNA.fai - > $output/bwa/mutation.readS.bam
        else
            echo "Unsupported type: {wildcards.TYPE}" && exit 1       
        fi

        ## Using bwa mapping result to remove false positive sites...
        if [[ $type_end == 'paired' ]]; then
            perl $Bin/bwa.fre.filt.pl $output/bwa/mutation.readP.bam $output/bwa/mutation.readS.bam $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.txt $annovar_humandb_refGene $output/readsvnum.txt $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.txt
        elif [[ $type_end == 'single' ]]; then
            perl $Bin/bwa.fre.filt_SE-1.pl $output/bwa/mutation.readS.bam $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.txt $annovar_humandb_refGene $output/readsvnum.txt $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.txt
        else
            echo "Unsupported type: {wildcards.TYPE}" && exit 1       
        fi
        """



rule PS81_2__Qiu2016_filter_variants____part3__sort_and_annovar_and_splicing_and_basicstat:
    input:
        flag_PS81_2__part2="result/PS81_2__Qiu2016_filter_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{SNPDB}/finished.part2__bwa_realignment_and_filtering"
    output:
        flag_PS81_2__part3=touch("result/PS81_2__Qiu2016_filter_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{SNPDB}/finished.part3__sort_and_annovar_and_splicing_and_basicstat")
    params:
        bam_directory="result/PS15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        call_directory="result/PS81_1__Qiu2016_call_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/",
        filter_directory="result/PS81_2__Qiu2016_filter_variants/{TYPE}/{DATASET}/{SAMPLE}/__merged__/{FASTP_PARAMETERS}/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{INDEXER}/{INDEXER_PARAMETERS}/{ALIGNER}/{ALIGNER_PARAMETERS}/{VARIANT_CALLER}/{VARIANT_CALLER_PARAMETERS}/{dbSNP_VERSION}/{dbSNP_SUBSET}/{SNPDB}/",
        perl_scripts_directory="external/Qiu2016.rep/RNA_editing_pipeline/bin/",
        annovar_directory="tools/annovar/",
        hg19_alu_bed="external/UCSC.Table.Browser/hg19_alu.bed",
        hg38_alu_bed="external/UCSC.Table.Browser/hg38_alu.bed"
    threads:
        1
    shell:
        """
        Bin="{params.perl_scripts_directory}"
        output={params.filter_directory}
        annovar={params.annovar_directory}
        ## RNA editing sites sorting...
        sort -k1.4n -k2n $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.txt > $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.sort.txt
        mv $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.sort.txt $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.txt
        
        annovar_humandb_refGene=""
        alu=""
        if [[ {wildcards.CONTIGS_FASTA_FILENAME} == "hg19.fa" ]] ; then
            annovar_humandb_refGene={params.annovar_directory}/humandb_hg19/hg19_refGene.txt
            alu={params.hg19_alu_bed}
        elif [[ {wildcards.CONTIGS_FASTA_FILENAME} == "hg38.fa" ]] ; then
            annovar_humandb_refGene={params.annovar_directory}/humandb_hg38/hg38_refGene.txt 
            alu={params.hg38_alu_bed}
        else 
            echo "Unsupported CONTIGS_FASTA_FILENAME: {wildcards.CONTIGS_FASTA_FILENAME}" && exit 1
        fi
        

        ## Annovar annotation
        ## Candidate RNA editing sites annotation...
        perl $Bin/annovar.pl $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.txt $annovar $annovar_humandb_refGene $alu $output/RNA_editing.sites.annotation.txt

        ## Splicing filter;
        ## Splicing region filtering...
        perl $Bin/alu.splicing.filter.pl $output/RNA_editing.sites.annotation.txt $output/RNA_editing.sites.annotation.splicing.filter.txt

        ## Basic information statistics;
        ## Basic information counting...
        perl $Bin/stat.pl $output/RNA_editing.sites.annotation.splicing.filter.txt $output/TypeDist.stat.txt $output/GenomeDist.stat.txt
        """
