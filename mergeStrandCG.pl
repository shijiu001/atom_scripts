use warnings;
use strict;


#chr12	3000068	+	0	0	CG	CGG
while(<STDIN>)
{
	chomp;
	my @l1=split;
	my $l2=<STDIN>;
	my @l2=split("\t",$l2);
	if($l1[1] != $l2[1]-1)
	{die "coor wrong"}
	else
	{
		my  $c=$l1[3]+$l2[3];
		my  $t=$l1[4]+$l2[4];
		if($c+$t==0)
		{print $l1[0],"\t",$l1[1]-1,"\t",$l2[1],"\tNA\t",$c,"\t",$t,"\n";}
		else
		{
		print $l1[0],"\t",$l1[1]-1,"\t",$l2[1],"\t",$c/($c+$t),"\t",$c,"\t",$t,"\n";
		}
	}

}
