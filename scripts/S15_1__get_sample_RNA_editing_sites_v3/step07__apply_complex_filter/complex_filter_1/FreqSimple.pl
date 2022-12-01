use strict;

if (@ARGV != 3){
	die "need to provide 1 input: name\n";
}

my ($simple,$input,$output) = @ARGV;

my $bedtools = "intersectBed";


open (f1,"$input") || die "Error $input";
open (o1,">$input.temp") || die "Error";
while (<f1>)
{
	$_=~s/\s+$//;
	my @a = split /\t/,$_;
	if ($a[5] > 0.1){print o1 "$a[0]\t$a[1]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\n";}
}
close f1;
close o1;

my $command = "$bedtools -a $input.temp -b $simple -v > $input.temp2";
system ("$command");
open (f1,"$input.temp2") || die "Error";
open (o1,">$output") || die "error write";
while (<f1>)
{
	my @a = split /\s+/,$_;
	print o1 "$a[0]\t$a[1]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\n";
}
close f1;
close o1;
system "rm $input.temp";
system "rm $input.temp2";
