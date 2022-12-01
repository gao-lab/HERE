use strict;
use threads;
use threads::shared;
use File::Basename;

my $Reference=$ARGV[0];
my $file = $ARGV[1];
my $bam = $ARGV[2];
my $outfile = $ARGV[3];
my $process = $ARGV[4];

my $PerlScript = dirname($0)."/BLAT_candidates.pl";
my $outdir = "_blatFiles";
if (!-e "$file$outdir"){mkdir "$file$outdir";}

open (f1,"$file") || die "Error open $file\n";
my $totalsite;
while(<f1>){$totalsite++;}
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
		my $params = join "\t",$file,$bam,$j;
		threads->new(\&Blat,$params);
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

open (o1,">$outfile") || die "Error";
for (my $i = 1;$i <= $process;$i++)
{
	open (f1,"$file$outdir/homo.TMP$i") || die "Error";
	while (<f1>){print o1 "$_";}
	close f1;
}
close o1;

sub Blat()
{
	my $params = shift;
	my ($input,$recal,$i) = split /\t/,$params;
	print "TMP$i\n";
	open (f1,"$input") || die "Error $input";
	open (o1,">$file$outdir/homo.$i.vcf") || die "Error $file$outdir/homo.$i.vcf";
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
	system "perl $PerlScript $Reference $file$outdir/homo.$i.vcf $recal $file$outdir/homo.TMP$i";
}
