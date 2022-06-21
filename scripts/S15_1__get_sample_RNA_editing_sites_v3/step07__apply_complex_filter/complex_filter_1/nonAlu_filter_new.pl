#use warnings;
use strict;

if (@ARGV != 4) {
	die "need to provide 1 inputs: shared variant file and outputfile repAlu/nonRep name\n";
}
my ($repeat,$inputfile, $outputfile1,$outputfile2) = @ARGV;
my $bedtools = "intersectBed";


############### Step 1: transfer #################
open (f1,"$inputfile") || die "error $inputfile";
open (o1,">$inputfile.temp") || die "error";
while (<f1>)
{my @a = split /\s+/,$_;print o1 "$a[0]\t$a[1]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\n";}
close f1;
close o1;
############### Step 2 : bedtools ###################
my $command = "$bedtools -a $inputfile.temp -b $repeat -wa > $outputfile1.temp";
system ("$command");
my %hash;
open (f1,"$outputfile1.temp") || die "Error";
open (o1,">$outputfile1") || die "error";
while (<f1>)
{
	my @a = split /\s+/,$_;
	my $b = join "\t",$a[0],$a[1];
	if ($hash{$b} eq "")
	{
		print o1 "$a[0]\t$a[1]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\n";
		$hash{$b} = 1;
	}
}
close f1;
close o1;

$command = "$bedtools -a $inputfile.temp -b $repeat -v > $outputfile2.temp";
system ("$command");
my %hash2;
open (f1,"$outputfile2.temp") || die "Error";
open (o1,">$outputfile2") || die "error";
while (<f1>)
{
	my @a = split /\s+/,$_;
	my $b = join "\t",$a[0],$a[1];
	if ($hash2{$b} eq "")
	{
		print o1 "$a[0]\t$a[1]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\n";
		$hash2{$b} = 1;
	}
}
close f1;
close o1;
system "rm $inputfile.temp";
system "rm $outputfile1.temp";
system "rm $outputfile2.temp";
