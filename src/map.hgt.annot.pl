#ReadID  HChrom  HPos    HCiGar  Hseq    Virus   VPos    VCiGAR  VSeq
#HWI-ST590:243:C371BACXX:1:1101:1224:27485       4       180829216       51M     NGACCATCCTGGCCAACATGGTGAAACCCTGTCTCTACTAAAAATACAAAA     TorquetenomidivirusDNAisolate_AB303557  2631    32M19S  CCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAGAATGTTGAATATTNGCCC
#chr10   118491630       HSPA12A

open FH, "$ARGV[0]" or die "";
%calls;
while(my $l = <FH>)	{
	chomp $l;
	@a=split(/\t/,$l);
	$id="$a[0]:$a[1]";
	if ($a[3] == -1)	{	$a[3] = "NA";}
	if ($a[2] eq "\.")	{	$a[2] = "NA";}
	$calls{$id}="$a[2]\t$a[3]";
	
}
close FH;
open FH, "$ARGV[1]" or die "";
$head=<FH>;
chomp $head;
print "GeneName\tDistToGene\t$head\n";
while(my $l = <FH>)     {
        chomp $l;
        @a=split(/\t/,$l);
	$id="$a[0]:$a[1]";
	if (defined $calls{$id})	{
		print "$calls{$id}\t$l\n";
	}
	#else	{
	#	print "NA\tNA\t$l\n";
	#}
}
close FH;
		

