library("data.table")


fread("./report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/combined.RNA.DNA.comparison.dt.csv.gz") -> combined.RNA.DNA.comparison.dt
fread("./report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/all.RNA.alignments.snpEff.vcf.reformatted.dt.gz", select=c("CHROM", "POS", "REF", "ALT")) -> all.RNA.alignments.simple.dt

{

    copy(combined.RNA.DNA.comparison.dt) -> .;
    .[is.DNA.detected==FALSE] -> .;
    .[event.summary %in% c("A>G", "A>G;T>C")] -> .;
    .[RNA.edit.type == "filtered"] -> .;
    ## reduce to edits (collapse samples)
    .[, list(CHROM, POS)] -> .;
    unique(.) -> .;
    ## get genomic info
    merge(x=., y=all.RNA.alignments.simple.dt, by=c("CHROM", "POS"), all.x=TRUE, all.y=FALSE) -> .;
    ## convert to bed
    .[, list(chrom=CHROM, chromStart=POS-1, chromEnd=POS-1+1, name=paste(sep="", REF, "_", ALT), score=".")] -> .;
    .[, strand:=c("A_G"="+", "T_C"="-")[name]] -> .;
    ## 
    . -> all.filtered.RNA.edits.genomic.info.bed.dt
    
}

fwrite(all.filtered.RNA.edits.genomic.info.bed.dt, "./report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/all.filtered.RNA.edits.genomic.info.bed", row.names=FALSE, col.names=FALSE, sep="\t")

system("bedtools slop -i report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/all.filtered.RNA.edits.genomic.info.bed -g external/contigs/hg38.fa.fai -b 3 > report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/all.filtered.RNA.edits.genomic.info.b3.extended.bed")

system("bedtools getfasta -s -fi external/contigs/hg38.fa -bed report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/all.filtered.RNA.edits.genomic.info.b3.extended.bed -fo report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/all.filtered.RNA.edits.genomic.info.b3.extended.fasta")

system("grep -v '>' report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/all.filtered.RNA.edits.genomic.info.b3.extended.fasta | tr '[a-z]' '[A-Z]' > report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/all.filtered.RNA.edits.genomic.info.b3.extended.sequence.lines")

system("cp -r tools/tsl report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/")

date()
system("cd report.ver2/pipeline.validation/210203-GSE144296.A375.12.variants/tsl/cgi-bin && chmod u+x ./tsl && chmod u+x ../pvalue/pvalue && ./tsl -P ../../all.filtered.RNA.edits.genomic.info.b3.extended.sequence.lines -N ../../../../../result/A01_5__plot_motif/210215-sixth-dataset/201221-fifth-phenotype-collection/temp_result_directory/background.7mer.lines -K N -O ../../all.filtered.RNA.edits.genomic.info.b3.extended.sequence.lines.tsl.plot.png -x -y -C nucleo_weblogo -U cm -W 8 -H 4 -R 300  && cd -")
date()
