library("data.table")
library("magrittr")

single.contig.Alu.only.snpEff.annotated.annotated.records.only.vcf.gz.filename <- snakemake@input[["single_contig_Alu_only_snpEff_annotated_annotated_records_only_vcf_gz_filename"]]
if (FALSE){
    single.contig.Alu.only.snpEff.annotated.annotated.records.only.vcf.gz.filename <- "result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/chr10.contigs.Alu.only.snpEff.annotated.annotated.records.only.vcf.gz"
}

## Here some annotation might contain 'noncoding transcript exon'; this is because the editing site falls onto a protein-coding transcript and a noncoding transcript at the same time. The downstream sorting of sites will assign it to the correct category if we only consider protein-coding transcripts and do not consider splicing sites, which is exactly what this manuscript wants to do (note that 3'UTR, 5'UTR and |p. are all features exclusive to protein-coding transcripts; also note that 5'UTR premature start codon has no |p. and is assigned to the 5'UTR subset, which is reasonable; annotations with |p. are synonymous, missense, stop lost and their combination with other annotations).

single.contig.Alu.only.snpEff.annotated.annotated.records.only.vcf.reformatted.dt <- fread(cmd=paste(sep="", "zcat ", single.contig.Alu.only.snpEff.annotated.annotated.records.only.vcf.gz.filename, " | grep -v '^##'"), header=TRUE) %>%
    setnames("#CHROM", "CHROM") %>%
    {.[, Annotation.subset:="intron"]} %>%
    {
        cat(date(), "looping over terms...\n");
        for (term in c("splice", "3_prime_UTR", "5_prime_UTR", "\\|p\\.")){
            cat(date(), "  updating ", term, "\n");
            .[grepl(term, INFO), Annotation.subset:=term]
        }
        .[Annotation.subset=="\\|p\\.", Annotation.subset:="CDS"]
    }

for (term in c("intron", "3_prime_UTR", "5_prime_UTR", "CDS")) {
    fwrite(single.contig.Alu.only.snpEff.annotated.annotated.records.only.vcf.reformatted.dt[Annotation.subset==term, list(CHROM, POS-1, POS)], snakemake@output[[paste(sep="", "single_contig_Alu_only_snpEff_annotated_annotated_records_only_", term, "_bed_filename")]], row.names=FALSE, col.names=FALSE, sep="\t")
}
