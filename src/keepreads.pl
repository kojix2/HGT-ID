open BAM , "samtools view -h $ARGV[0] | " or die "";
#open OUTBAM, "  | samtools view -bS - >$ARGV[1]" or die ""; 
#open OUTBAM, " >$ARGV[1]" or die ""; 
### trying to remove the reads which are unpaired and bam falg is not set correclty
$prev="";
$prev_read="";
while(my $l = <BAM>)	{
	chomp $l;
	#if ($l =~ m/^@/) {print OUTBAM "$l\n"; next;}
	if ($l =~ m/^@/) {print  "$l\n"; next;}
	@a=split(/\t/,$l);
	@b=split(/\//,$a[0]);
	if ($b[0] eq $prev)	{
		#print OUTBAM "$prev_read\n";
		#print OUTBAM "$l\n"; 
		#print  "$prev_read\n";
		@a1=split(/\t/,$prev_read);
		@b1=split(/\//,$a1[0]);
		### to fix the wrong flag issue for rna seq samples
		### 8/25 liver TCGA sample
		if ($a1[1] == 4 )	{ $a1[1] = 69;}
		print  "$b1[0]\t". join ("\t",@a1[1..10]) . "\n";
		if ($a[1] == 4 )	{ $a[1] = 133;}
		print  "$b[0]\t". join ("\t",@a[1..10]) . "\n";
		#print "no\n";
	}	else	{
		#print "yes\n";
		$prev_read=$l;
		$prev=$b[0];
	}
}
close BAM;
close OUTBAM;	
		