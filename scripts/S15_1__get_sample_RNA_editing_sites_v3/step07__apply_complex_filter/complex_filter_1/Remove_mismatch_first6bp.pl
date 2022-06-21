#script that removes mismatches in first 6 bp of reads that are caused by imperfect priming while making cDNA from RNA

use warnings;
use strict;

if (@ARGV != 3) {
	die "need to provide 3 input:Variant list, INDEXED BAM alignment file and output file name\n";
}
my ($inputfile, $bamfile, $outputfile) = ($ARGV[0], $ARGV[1], $ARGV[2]);

my $minbasequal = 25;
my $minmismatch = 1;
my $samtools="samtools";

open (my $INPUT , "<", $inputfile) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $outputfile);

while (<$INPUT>) {
	chomp;
	my @fields = split;
	my $TEMPNAME = join '', $outputfile,'_tmp';
	my ($chr, $position) = ($fields[0], $fields[1]);
	my $bamposition = join '', $chr,':',$position,'-',$position;
	system("$samtools view $bamfile $bamposition > $TEMPNAME"); #get all reads overlapping position
	my $editnuc = $fields[4];
	my ($newcov, $newmismatch) = (0,0);

	open(my $TMPFILE, "<", $TEMPNAME);
	while (<$TMPFILE>) { #loop through the reads
		chomp;
		my @bamfields = split;
		my ($alignment, $readstart, $cigar, $sequence, $qualities) = ($bamfields[1], $bamfields[3], $bamfields[5], $bamfields[9], $bamfields[10]);
		my @sequencebases = split(//,$sequence);
		my @qualscores = split(//,$qualities);
		my ($currentpos, $readpos) = ($readstart,1);
		my $base_readpos;
		my @cigarnums = split(/[MIDNSHP]/, $cigar);
		my @cigarletters = split(/[0-9]+/,$cigar);
		shift @cigarletters;

		for (my $i = 0; $i < @cigarnums; $i++) {
			if ($cigarletters[$i] =~ m/[ISH]/) {
				$readpos = $readpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/[DN]/) {
				$currentpos = $currentpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/M/) {
				for (my $j = 0; $j < $cigarnums[$i]; $j++) {
					$base_readpos = $readpos if ($currentpos == $position); #find position in read that overlaps the candidate site
					$currentpos++;
					$readpos++;
				}	
			}
		}
		next unless $base_readpos;
		my $revstrand = 0;
		$revstrand = 1 if ($alignment & 16);
		if (($revstrand == 0 && $base_readpos > 6) || ($revstrand == 1 && $base_readpos < $readpos - 5)) { #ensure that this position is not in first 6 bp of read
			if (ord($qualscores[$base_readpos-1]) >= $minbasequal+33) {  #SANGER format!
				$newcov++;
				$newmismatch++ if ($sequencebases[$base_readpos-1] eq $editnuc); #check if this read is a "wild type" or "mismatched" read	
			}
		}	
	}

	close($TMPFILE);
	#system("rm $TEMPNAME");

	if ($newmismatch >= $minmismatch) { #only output candidates that have mismatched reads remaining
		my $varfreq = sprintf("%.3f", $newmismatch/$newcov);
		print $OUTPUT "$fields[0]\t$fields[1]\t$newcov,$newmismatch\t$fields[3]\t$fields[4]\t$varfreq";
		for (my $i = 6; $i < @fields; $i++) {
			print $OUTPUT "\t$fields[$i]";
		}
		print $OUTPUT "\n";
	}
}
close $INPUT;	
close $OUTPUT;
