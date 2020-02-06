open FH, "$ARGV[0]" or die "";

open FA, "$ARGV[1]" or die "";

%contig;

while(my $l = <FA>)	{
	@a=split(/\t/,$l);
	$contig{$a[0]}=1;
#	print "$a[0]\n";	
}
close FA;

while(my $l = <FH>)	{
	chomp $l;
	@a=split(/\t/,$l);
		if(defined($contig{$a[2]}) && defined($contig{$a[6]}))	{
		}
		else	{
			print "$l\n";
		}	 	
}
close FH;		
