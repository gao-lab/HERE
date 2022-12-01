
#script used to filter out variants called due to intronic alignment of RNA-SEQ reads
use warnings;
use strict;

if (@ARGV != 3) {
	die "need to provide 3 inputs:list of variants gene annotation file in UCSC format and output file name\n";
}
my ($inputfile, $genefile , $output) = @ARGV;
my %genehash;

open (GENEFILE, $genefile) or die "error opening gene annotation file: $!\n";
open (OUTPUT, ">$output") or die "error opening output file: $!\n";
while (<GENEFILE>) {
	chomp;
	my @fields = split;
	my $chr = $fields[2];
	push(@{$genehash{$chr}},$_);
}
close GENEFILE;

open (INPUTFILE, $inputfile) or die "error opening input variant list: $!\n"; #open file of candidate editing sites
while (<INPUTFILE>) {
	chomp;
	my $line = $_;
	my @fields = split(/\t/);
	my ($chrom, $position, $found, $chromfound, $splice) = ($fields[0], $fields[1], 0, 0, 0);
	my ($gene, $strand, $gene2) = ('','','');
	foreach my $geneline (@{$genehash{$chrom}}) { #for each variant, loop through gene annotation file and check if variant is in intronic sequence near an intron/exon boundary
		my @fieldsref = split(/\t/,$geneline);
		my ($chromref, $txstart, $txend) = ($fieldsref[2], $fieldsref[4], $fieldsref[5]);
		$chromfound = 1;
		if (int($txstart) <= int($position) && int($txend) >= int($position)) {
			$found = 1;
			$strand = join '',$strand ,',' ,$fieldsref[3];
			$gene = join '',$gene ,',' ,$fieldsref[1];
			$gene2 = join '',$gene2 ,',' ,$fieldsref[12];
			my $exoncount = int($fieldsref[8]);
			my @exonstarts = split(/,/, $fieldsref[9]);
			my @exonends = split(/,/, $fieldsref[10]);
			for (my $i = 0; $i < $exoncount; $i++) { #for each exon, check if variant lies in it
				$splice = 1 if ((int($exonstarts[$i])-4 < int($position) and int($exonstarts[$i])+1 > $position) or (int($exonends[$i]) < $position and int($exonends[$i])+4 >= $position)); #make sure variant is not within 4 bp to intronic side of ANY intron/exon boundary		
			}
		} 
	}
	if ($splice == 0) { #if variant is not intron near splice boundary, print it out
		print OUTPUT "$line\n"; 
	} 
}
close INPUTFILE;
