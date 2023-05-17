use warnings;
use strict;
my @args=@ARGV;
#my $cg_list=$args[-1];
my %c;my %t;
for(my $i=0;$i<=$#args;$i++)
{
	open IN,$args[$i];
	while(<IN>) 
	{
		chomp;
		my @l=split;
		my $c=$l[-2];
		my $t=$l[-1];
		my $pos=$l[1];
		$c{$l[0]."\t".$pos}+=$c;
		$t{$l[0]."\t".$pos}+=$t;
		
	}
	close IN;
}

open IN1,$args[0];
while(<IN1>)
{
	chomp;
	my @l=split;
	my $b=0;
	my $ml="NA";
	if(exists $c{$l[0]."\t".$l[1]} || exists $t{$l[0]."\t".$l[1]})
	{
		$b=$c{$l[0]."\t".$l[1]}+$t{$l[0]."\t".$l[1]};
		$ml=$c{$l[0]."\t".$l[1]}/$b if $b != 0;
		if($b != 0)
		{
		print "$l[0]\t",$l[1],"\t$l[2]\t",sprintf("%.4f", $ml),"\t",$c{$l[0]."\t".$l[1]},"\t",$t{$l[0]."\t".$l[1]},"\n";
		}
	}
}
close IN1;
