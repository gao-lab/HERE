use strict;
use threads;
use threads::shared;
use File::Basename;

my $Reference=$ARGV[0];
my $input = $ARGV[1];
my $geneAnno = $ARGV[2];
my $process = $ARGV[3];

my $outdir = "$input\_rmSJandHome/";
if (!-e $outdir){mkdir $outdir;}

## my $transcript = "/lustre4/gaog_pkuhpc/users/dingy/Project.collaborated/Project.embryo.I/200827-results/200827-lihao-scripts-modified/RNAediting_scripts_ouyang/Transcript/";
my $HomoFilter = dirname($0)."/RemoveHomoNucleotides.pl";
my $SJFilter = dirname($0)."/Filter_intron_near_splicejuncts.pl";

open (f1,"$input.sim.vcf") || die "Error";
my $totalsite;
while (<f1>){$totalsite++;}
close f1;
my $eachNum = int($totalsite/$process)+1;

my $j = 0;
my $thread;
while (1)
{
	last if ($j >= $process);
	while (scalar(threads->list()) < $process)
	{
		$j++;
		my $params = join "\t",$input,$geneAnno,$j;
		threads->new(\&RmCommand,$params);
	}
	foreach $thread(threads->list(threads::all))
	{
		if ($thread->is_joinable())
		{
			$thread->join();
		}
	}
}

foreach $thread(threads->list(threads::all)){$thread->join();}

open (o1,">$input.rmSJandHomo.txt") || die "Error";
for (my $i = 1;$i <= $process;$i++)
{
	open (f1,"$outdir$i.filter.Homo") || die "Error";
	while (<f1>){print o1 "$_";}
	close f1;
}
close o1;

sub RmCommand()
{
	my $params = shift;
	my ($input,$geneAnno,$i) = split /\t/,$params;
	print "TMP$i\n";
	open (f1,"$input.sim.vcf") || die "Error";
	open (o1,">$outdir$i.sim.temp") || die "$outdir$i.sim.temp";
	my $num = 0;
	my $n;
	while (<f1>)
	{
		if ($num < ($i)*$eachNum && $num >= ($i-1)*$eachNum)
		{$n++;print o1 "$_";}
		$num++;
	}
	close f1;
	close o1;
	#filter of Splice Junction
	system ("perl $SJFilter $outdir$i.sim.temp $geneAnno $outdir$i.filter.SpliceJunction");
	#filter of homopolymer runs of 5nt
	system ("perl $HomoFilter $Reference $outdir$i.filter.SpliceJunction $outdir$i.filter.Homo");
}
