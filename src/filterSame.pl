#Hchr	Hstart	Hend	Virus	Vstart	VStop	Reads
#chrGL000214.1	21860	21860	HepatitisBvirus_21326584	1550	1550	11
#chrGL000214.1	33522	33535	chr1	142621855	142622135	5

open FH, "$ARGV[0]" or die "";

open FA, "$ARGV[1]" or die "";

%contig;

while(my $l = <FA>)	{
	@a=split(/\t/,$l);
	$contig{$a[0]}=1;
}
close FA;

while(my $l = <FH>)	{
	chomp $l;
	if ($. == 1)	{print "$l\n";}
	else	{
		@a=split(/\t/,$l);
		if(defined($contig{$a[0]}) && defined($contig{$a[3]}))	{
		}
		else	{
			print "$l\n";
		}	 	
	}
}
close FH;


