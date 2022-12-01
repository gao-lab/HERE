#BLAT mismatched reads to ensure unique mapping!

#use warnings;
use strict;

if (@ARGV != 4) {
	die "need to provide 3 input:Variant list, INDEXED BAM alignment file and output file name\n";
}
my ($Reference,$inputfile, $bamfile, $outputfile) = @ARGV;

my $minbasequal = 25;
my $minmismatch = 1;
my $scorelimit = 0.95; #second blat hit score must be less than 95% of first hit!


open (my $INPUT , "<", $inputfile) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $outputfile);

my $fafile = join '', $inputfile, '.fatmp';
my $pslfile = join '', $inputfile, '.psltmp';
open (FAFILE,">$fafile");

while (<$INPUT>) {
	chomp;
	#if($_ !~ /^$ARGV[3]\t/){next;}
	my $inputline = $_;
	my @fields = split;
	my $TEMPNAME = join '', $outputfile,'_tmp';
	my ($chr, $position) = ($fields[0], $fields[1]);
	my $bamposition = join '', $chr,':',$position,'-',$position;
	system("samtools view $bamfile $bamposition > $TEMPNAME"); #get reads overlapping candidate site
	my $editnuc = $fields[4];
	my $newmismatch = 0;
	my $mismatchreadcount = 0;

	open(my $TMPFILE, "<", $TEMPNAME);
	while (<$TMPFILE>) {
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
				for (my $j = 0; $j < $cigarnums[$i]; $j++) { #for each read, determine if it contains the mismatch at the candidate site
					$base_readpos = 1 if ($currentpos == $position && $sequencebases[$readpos-1] eq $editnuc && ord($qualscores[$readpos-1]) >= $minbasequal+33); #SANGER format!
					$currentpos++;
					$readpos++;	
				}	
			}
		}
		if ($base_readpos) {
			print FAFILE ">$chr\-$position\-$mismatchreadcount\n$sequence\n"; #make fasta file of all mismatched reads
			$mismatchreadcount++;
		}
	}
	close $TMPFILE;
	system("rm $TEMPNAME");
}
close $INPUT;

system("blat -stepSize=20 -repMatch=2253 -minScore=20 -minIdentity=0 -noHead $Reference $fafile $pslfile");
system("rm $fafile");



open(PSL, "<", $pslfile );
my %pslhash;
while(<PSL>) { #parse PSL outputfile: for each mismatched read that was blatted, hash together all the blat hits 
	chomp;
	my @pslfields = split;
	my $name = $pslfields[9];
	my $blatscore = join '@',$pslfields[0],$pslfields[13],$pslfields[17],$pslfields[18],$pslfields[20];
	if ($pslhash{$name}) {
		$pslhash{$name} = join '-', $pslhash{$name}, $blatscore;
	} elsif ( !$pslhash{$name}) {
		$pslhash{$name} = $blatscore; 
	}
}
close PSL;	

my %sitehash;
my %discardhash;
foreach my $pslkey (keys %pslhash) { #look at each blatted read one by one
	my @splitkey = split(/\-/, $pslkey);
	my $site = join '_', $splitkey[0],$splitkey[1]; #assign this read to its corresponding candidate site
	my @psllines = split(/\-/, $pslhash{$pslkey}); #separate out the blat hits for each read
	my $largestscore = 0;
	my $largestscoreline = $psllines[0];
	my @scorearray;
	foreach my $scoreline (@psllines) { #find the largest scored blat hit in the genome
		my @scoresarray = split(/\@/,$scoreline);
		my $linescore = $scoresarray[0];
		push(@scorearray,$linescore);
		if ($linescore > $largestscore) {
			$largestscoreline = $scoreline;
			$largestscore = $linescore;
		}
	}
	@scorearray = sort {$b <=> $a} @scorearray; #sort the blat hits by score
	$scorearray[1] = 0 unless ($scorearray[1]);
	my @splitlargestline = split(/\@/,$largestscoreline);
	my $overlapfound = 0;
	if ($splitlargestline[1] eq $splitkey[0] && $scorearray[1] < ($scorearray[0]*$scorelimit)) { #ensure that score of second blat hit is less than 95% of first hit
		my ($numblocks, $blocksizes, $blockstarts) = ($splitlargestline[2],$splitlargestline[3],$splitlargestline[4]);
		my @blocksizes = split(/\,/,$blocksizes);
		my @blockstarts = split(/\,/,$blockstarts);
		for (my $i = 0; $i < $numblocks; $i++) {
			my $startpos = $blockstarts[$i]+1;
			my $endpos = $blockstarts[$i] + $blocksizes[$i];
			$overlapfound = 1 if ($splitkey[1] >= $startpos && $splitkey[1] <= $endpos); #ensure that first blat hit overlaps the candidate editing site
		}
		if ($sitehash{$site} && $overlapfound) { #check if this read passes the blat criteria
			$sitehash{$site}++;
		} elsif (!$sitehash{$site} && $overlapfound) {
			$sitehash{$site} = 1;
		}
	}
	unless ($overlapfound) { #check if this read has failed the blat criteria
		if ($discardhash{$site}) {
			$discardhash{$site}++;
		} elsif (!$discardhash{$site}) {
			$discardhash{$site} = 1;
		}
	}
}

open (SITES2, "<", $inputfile ) or die "error opening inputfile: $!\n"; #open input file again and check if each candidate passes the blat filtering
while(<SITES2>) {
	chomp;
	my @fields = split;
	my $inputline = $_;
	my ($cov,$oldalter) = split(/\,/,$fields[2]);
	my $name = join '', $fields[0],'_',$fields[1];
	if ($sitehash{$name}) {
		my $newalter = $sitehash{$name};
		my $discardnum;
		if ($discardhash{$name}) {
			$discardnum = $discardhash{$name};
		} else {
			$discardnum = 0;
		}
		my $newcov = $cov - ($oldalter - $newalter);
		my $neweditfreq = sprintf("%.3f", $newalter/$newcov);
		print $OUTPUT "$fields[0]\t$fields[1]\t$newcov,$newalter\t$fields[3]\t$fields[4]\t$neweditfreq\n" if ($newalter >= $minmismatch && $newalter > $discardnum); #make sure that number of read passing blat criteria is greater than number that fail
	}
}
close SITES2;
close $OUTPUT;
#system("rm $pslfile");
