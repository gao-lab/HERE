#use warnings;
use strict;

if (@ARGV != 4) {
	die "need to provide 4 input:Variant list, reference variant file and output file name\n";
}
my ($cdna,$inputfile, $refvar, $outputfile) = @ARGV;

open (INPUT , "<", $inputfile) or die "error opening inputfile: $!\n";
open (REF , "<", $refvar) or die "error opening ref: $!\n";
open (CDNA , "<", $cdna) or die "error opening cdna: $!\n";
open (OUTPUT, ">", $outputfile);

my @a;
my %temp;
my %cd;
my $t;

while (@a = split /\t/,<CDNA>){
	if ($a[0] !~ /#/){$cd{$a[3]} = 1;}
}
close CDNA;

while (@a = split /\t/,<REF>)
{
	if ($a[0] !~ /#/){
		my $c = join "\t",$a[0],$a[1];
		if ($cd{$a[2]} == 1 || $a[2] !~ /rs/){$temp{$c} = 1;} 
	}
}
close REF;

while (<INPUT>){
	@a = split /\t/,$_;	
	$t = join "\t",$a[0],$a[1];
	if ($temp{$t} == 1){
		print OUTPUT "$_";
	}
}
close INPUT;
close OUTPUT;
