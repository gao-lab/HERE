#remove variants that lie in homopolymer repeat (analogous to GATK "HRun" filter)

use warnings;
use strict;

if (@ARGV != 3) {
	die "need to provide 2 input:editing site file and output file name\n";
}
my ($Reference,$inputfile, $outputfile) = @ARGV;
my ($leftbuffer, $rightbuffer) = (4,4); #sequence buffer length on both sides

my $tmpfile = join '', $inputfile, '.nuctmp';
my $tmpfafile = join '', $inputfile, '.tmpnuctmp';


open (SITES, "<", $inputfile ) or die "error opening inputfile: $!\n"; #Read INPUTFILE
open (OUTPUT, ">", $outputfile);


while(<SITES>) { #Read INPUTFILE
	#if($_ !~ /^$ARGV[2]\t/){next;}
	print "$_";
	chomp;
	my $line = $_;
	my @fields = split;
	my ($chromname, $position, $editnuc) = ($fields[0], $fields[1], $fields[4]);
	my $homopol = 0;
	
	open (TMP, ">", $tmpfile);
	my ($startpos, $endpos) = ($position-$leftbuffer-1,$position+$rightbuffer);
	print TMP "$chromname\t$startpos\t$endpos\n";
	print "$chromname\t$startpos\t$endpos\n";
	close TMP;

	system("fastaFromBed -fi $Reference -bed $tmpfile -fo $tmpfafile"); #get sequence surrounding candidate position
	my $seq = '';
	open (TMPFAFILE, "<", $tmpfafile);
	while(<TMPFAFILE>) {
		chomp;
		next if ($_ =~ m/\>/);
		$seq = join '', $seq, $_; #get sequence flanking variant
	} 
	close TMPFAFILE;
	system("rm $tmpfafile");
	system("rm $tmpfile");
	my @splitseq = split(//,$seq);
	$editnuc = uc $editnuc;
	for (my $i = 0; $i < @splitseq; $i++) {
		$splitseq[$i] = uc $splitseq[$i];
	}
	$homopol = 1 if ($splitseq[0] eq $editnuc && $splitseq[1] eq $editnuc && $splitseq[2] eq $editnuc && $splitseq[3] eq $editnuc);
	$homopol = 1 if ($splitseq[1] eq $editnuc && $splitseq[2] eq $editnuc && $splitseq[3] eq $editnuc && $splitseq[5] eq $editnuc);
	$homopol = 1 if ($splitseq[2] eq $editnuc && $splitseq[3] eq $editnuc && $splitseq[5] eq $editnuc && $splitseq[6] eq $editnuc);
	$homopol = 1 if ($splitseq[3] eq $editnuc && $splitseq[5] eq $editnuc && $splitseq[6] eq $editnuc && $splitseq[7] eq $editnuc);
	$homopol = 1 if ($splitseq[5] eq $editnuc && $splitseq[6] eq $editnuc && $splitseq[7] eq $editnuc && $splitseq[8] eq $editnuc);
	print OUTPUT "$line\n" unless ($homopol); 
}

close SITES;
close OUTPUT;
