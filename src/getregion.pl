
use List::Util qw( min max );

# $scan=500;

@hold=();
$prev_chr="";
$prev_vchr="";
$prev_pos=0;
open FH, "$ARGV[0]" or die "";
$scan=$ARGV[1];
print "#Hchr\tHstart\tHend\tVirus\tVstart\tVStop\tReads\n";
while(my $l = <FH>)	{
	chomp $l;
	@a=split('\s+',$l);
	$cmp=$prev_pos + $scan;
	#print "$cmp\n";
	#print "$a[0]\t$a[1]\n";
	if ($a[0] eq $prev_chr && $a[1] <= $cmp && $a[2] eq $prev_vchr)	{
		#print "if\n";
		$prev_chr=$a[0];
		$prev_pos=$a[1];
		$prev_vchr=$a[2];
		$add=join(":",@a);
		push (@hold,$add);
		#print join ("\t",@hold) . "\n" ;	
		#<STDIN>;
	}
	else	{
		if (scalar(@hold) >= 1)	
		{
			#print "else\n";
			@b=split(/:/,$hold[0]);
			$hstart=$b[1];
			$hchr=$b[0];	
			$vchr=$b[2];
			@b=split(/:/,$hold[$#hold]);
			$hstop=$b[1];
			@reads=();
			$depth=0;	
			for($i=0;$i<=$#hold;$i++)	{
				@a1=split(/:/,$hold[$i]);
				push(@reads,$a1[$#a1 -1]);	
				$depth+=$a1[$#a1];
			}
			$min = min @reads;
			$max = max @reads;
			#print "$min\t$max\n";
			print "$hchr\t$hstart\t$hstop\t$vchr\t$min\t$max\t$depth\n";	
		}	
		@hold=();	
		$prev_chr=$a[0];
		$prev_pos=$a[1];
		$prev_vchr=$a[2];		
		$add=join(":",@a);
		push (@hold,$add);
		#print join ("\t",@hold) . "\n" ;	
	}	
}
if (@hold)	{
	@b=split(/:/,$hold[0]);
	$hstart=$b[1];
	$hchr=$b[0];
	$vchr=$b[2];	
	@b=split(/:/,$hold[$#hold]);
	$hstop=$b[1];
	@reads=();
	$depth=0;	
	for($i=0;$i<=$#hold;$i++)	{
		@a1=split(/:/,$hold[$i]);
		push(@reads,$a1[$#a1 -1]);	
		$depth+=$a1[$#a1];
	}
	$min = min @reads;
	$max = max @reads;	
	print "$hchr\t$hstart\t$hstop\t$vchr\t$min\t$max\t$depth\n";
}			
close FH;
