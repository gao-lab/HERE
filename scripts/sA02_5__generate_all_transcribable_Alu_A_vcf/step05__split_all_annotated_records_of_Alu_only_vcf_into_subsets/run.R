library("data.table")

single.contig.Alu.only.combined.snpEff.annotated.annotated.records.only.vcf.gz.filename <- snakemake@input[["all_contigs_Alu_only_snpEff_annotated_annotated_records_only_combined_vcf_gz_filename"]]
if (FALSE){
    all.contigs.Alu.only.combined.snpEff.annotated.annotated.records.only.vcf.gz.filename <- "result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/all.contigs.Alu.only.snpEff.annotated.annotated.records.only.combined.vcf.gz"
}


all.contigs.Alu.only.combined.snpEff.annotated.annotated.records.only.vcf.reformatted.dt <- fread(cmd=paste(sep="", "zcat ", all.contigs.Alu.only.combined.snpEff.annotated.annotated.records.only.vcf.gz.filename, " | grep -v '^##'"), header=TRUE, drop=c("QUAL", "FILTER")) %>%
    setnames("#CHROM", "CHROM") %>%
    {.[, Annotation.subset:="intron"]} %>%
    {
        cat(date(), "looping over terms...\n");
        for (term in c("splice", "3_prime_UTR", "5_prime_UTR", "\\|p\\.")){
            cat(date(), "  updating ", term, "\n");
            .[grepl(term, INFO), Annotation.subset:=term]
        }
        .[Annotation.subset:="\\|p\\.", Annotation.subset:="CDS"]
    }
