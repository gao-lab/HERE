#!/usr/bin/perl

use lib qw(BioPerl-1.6.923); # NEEDS BIOPERL
$| = 1;
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Spec;
use File::Basename;
use Bio::DB::Fasta;

my $length;
my $genomefile;
my $bedfile;
my $outname;
if(@ARGV)
{
    $genomefile=$ARGV[0];
    $bedfile=$ARGV[1];
    $length=$ARGV[2];
    $outname=$ARGV[3];


}

## my $outname = join "",$length,"bp";
open o1,"> $outname" || die "error";#">/home/ouyang/project3-RNAediting/RNAediting-protocol/SplicingJunction_Index/splice_junctions_$outname.txt"

my $INDIR = "./"; ##"/home/ouyang/data/hg19_chrs/";## "/home/ouyang/project3-RNAediting/RNAediting-protocol/SplicingJunction_Index/mm10_chrs/";
#directory of all human chromosomes as individual Fasta fileshg19-from-Renchao/hg19_chrs/

my @bedFiles = $bedfile;

my $WINDOWSIZE = $length;
my $GENOMEFILE = $genomefile;

my $JUNCTIONS;
foreach(@bedFiles){
	print STDERR $_."\n";
	run($INDIR.$_);
}

foreach my $key(sort(keys %$JUNCTIONS)){
	print o1 ">$key\n";
	print o1 $JUNCTIONS->{$key}."\n";
}



sub run{
	my $annotfile = shift;
	my $junctions;
	
	############
	# read genome file
	my $seqs;
	print STDERR "reading genome file...\n";
	my $seqDb =  Bio::DB::Fasta->new($GENOMEFILE);
	
	############
	# read annotation file
	my $annot = readAnnotFile($annotfile);
	#exit;
	
	############
        
	# go through all spliceforms
	foreach my $key(sort(keys %$annot)){
		
		my @exonStarts = @{$annot->{$key}->{exonStarts}};
		my @exonEnds = @{$annot->{$key}->{exonEnds}};
		my ($chromosome) = $key =~ /(.*):.*/;
		
		next if($#exonStarts == 0);

		for(my $i = 1; $i <= $#exonStarts; $i++){
			
			my ($lSeq,$lSeqPos) = getLseq(\@exonStarts,\@exonEnds, $i-1,$seqDb,$chromosome);
			my ($rSeq,$rSeqPos) = getRseq(\@exonStarts,\@exonEnds, $i,$seqDb,$chromosome);

			my $junctionId =  $chromosome."-".$lSeqPos."-".$rSeqPos;

			## print STDERR "$junctionId\n";
			$JUNCTIONS->{$junctionId} = $lSeq.$rSeq;
		}
		
	}
	
	#print "\n";
}

sub getLseq{
	my ($exonStarts, $exonEnds, $i,$seqDb,$chromosome) = @_;
	
	my $remaining = $WINDOWSIZE;
	my $lSeq = "";
	my @lSeqPosArr = ();
	
	while($remaining>0){
		if($i>=0){
			
			my $winStart = $exonEnds->[$i]-$remaining < $exonStarts->[$i] ? $exonStarts->[$i] : $exonEnds->[$i]-$remaining;
			#print "$winStart\n";
			$lSeq = uc($seqDb->get_Seq_by_id($chromosome)->subseq($winStart+1=>$exonEnds->[$i])).$lSeq;
			unshift @lSeqPosArr, ($winStart+1)."-".$exonEnds->[$i];
			$remaining -= ($exonEnds->[$i] - $winStart);
			$i -= 1;
		}
		else{
			return ($lSeq,join("-",@lSeqPosArr));
		}
	}
	return ($lSeq,join("-",@lSeqPosArr));	
}

sub getRseq{
	my ($exonStarts, $exonEnds, $i,$seqDb,$chromosome) = @_;
	
	my $remaining = $WINDOWSIZE;
	my $rSeq = "";
	my @rSeqPosArr = ();
	
	while($remaining>0){
		if($i<=scalar(@{$exonStarts})-1){
			my $winEnd = $exonStarts->[$i]+$remaining > $exonEnds->[$i] ? $exonEnds->[$i] : $exonStarts->[$i]+$remaining;
			$rSeq = $rSeq.uc($seqDb->get_Seq_by_id($chromosome)->subseq($exonStarts->[$i]+1=>$winEnd));
			push @rSeqPosArr, ($exonStarts->[$i]+1)."-".$winEnd;
			$remaining -= ($winEnd - $exonStarts->[$i]);
			$i += 1;
		}
		else{
			return ($rSeq,join("-",@rSeqPosArr));			
		}
	}
	return ($rSeq,join("-",@rSeqPosArr));		
}




sub readAnnotFile{
	my $annotfile = shift;
	my $annot;
	print STDERR $annotfile."\n";
	open(IN, "<".$annotfile);
	my @filecnt = <IN>;
	shift @filecnt; 
	foreach(@filecnt){
		my $line = $_;
		my @cnt = split("\t", $line);
		if($cnt[2] !~ /_/ )#filter the random chromosomes ,v1==/chr.*_.*/,note that the name of chromosome is "chr1" or just number "1"
		{
			$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5].$cnt[1]}->{exonStarts} = [split(",",$cnt[9])];
			$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5].$cnt[1]}->{exonEnds} = [split(",",$cnt[10])];
			$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5].$cnt[1]}->{gene} = $cnt[1];
		}
	}
	
	return $annot;
}


