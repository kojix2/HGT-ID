#ERR092850.522	83	chr12	5423805	70	3S86M1S	=	5423805	86	ATGAGGAAAAGGCTTTGCACATACTTAATGAGCTCTGAGTCAGTCTCATAAGTGGAATGTTCTAAGGAACCACCAGACAGGTAAAATAAN	D=DDDD=@DD?????:>>;:??;??=DD9DBDD@@DB@DDBD==DDBDDDDDDDD@@=@@@@@DD@@@@@@D@D@DBDDD<22/<.622#	LB:Z:177T	AS:i:56	UQ:i:56	NM:i:0	MD:Z:86	PQ:i:395	SM:i:70	AM:i:8	RG:Z:177T	PG:Z:novoalign
#SRR1609141.77670609	133	chr1	10001	0	*	=	10001	0	TTAGGGTTAGGGTTAGGGTTAAGGGTTGGGGTTAGGGTTAGGGTTAGGGTTAGGGTTGGGGTAAGGGTTAGGGTATGGGTTAGGTATGGG	a_[a\Y`Z`b]ecfdgfee_cae`e^HHYcY_HHWHaO^^Z_eM\\Q\W\bBB
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
#SRR1609151.1171895	133	chr1	10001	0	*	=	10001	0	TAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGGTAGGGTTAGGGTTAGGGTTAGGGGTAGGGGTAGGGGTAGGGGTAGGGT	a_beecccggfeeghc_`eeghiXbdddbbbfd`f^Wa_]eHV_egfMWZ\__
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

$count=0;
while(<>)	{
	chomp ;
	@_ = split(/\t/, $_);
	### length of the read
	$readlen=length(@_[9]);
	$cmp= $readlen/2;
	my @CIGAR = split(/([0-9]+[SMIDNHXP])/, @_[5]);
	#print "$CIGAR[$#CIGAR]\n"; 
	#map { push(@{$hash->{$2}}, $1) if (/(\d+)([SMIDNHXP])/) } @CIGAR;
	### first part
	if ($CIGAR[1] =~ m/S/)	{
		$tag=$CIGAR[1];
		$tag =~ m/(\d+)([SMIDNHXP])/;
		if ($1 >= $cmp)	{
			$len=scalar(@CIGAR);
			#print "length=$len\n";
			$base=substr(@_[9],1,$1);
			$qual=substr(@_[10],1,$1);	
			print "@_[0]_$count\t69\t@_[2]\t@_[3]\t0\t*\t=\t@_[3]\t0\t$base\t$qual\n";
			$base=substr(@_[9],$1);
                        $qual=substr(@_[10],$1);      
			#print "@CIGAR\n";
			#$Ci=join("a",@CIGAR[2..$CIGAR[$#CIGAR]]);
			$Ci=@_[5];
			$Ci =~ s/$CIGAR[1]//;
			if (@_[1] & 0x10){
				$bitflag=153;
			}else	{
				$bitflag=137;
			}
			print "@_[0]_$count\t$bitflag\t@_[2]\t@_[3]\t@_[4]\t$Ci\t@_[6]\t@_[3]\t0\t$base\t$qual\n";
		}
		$count++;
	}
#ERR092850.749	147	chr2	17933664	70	41M49S	=	17933601	-104	GAAATATTTTCAGTCCCTTGATTGAGTCCACGGATACAGAGTAAGACCCTGTCTCAAAAAAAAAAAAGGTTATCAATGAATTATTTGCCA	=AC8;+DDDD@CC3FFEFEFBGBCFFDGGCFGFEF==@??DAADA.DBDDFC>BAEHGHHHHHHHHHHHHHHHHHHHHEHFHHHHHHHHH	LB:Z:177T	AS:i:285	UQ:i:285	NM:i:0	MD:Z:41	PQ:i:320	SM:i:9	AM:i:9	RG:Z:177T	PG:Z:novoalign

	if ($CIGAR[$#CIGAR] =~ /S/)	{
		$tag=$CIGAR[$#CIGAR];
		$tag =~ m/(\d+)([SMIDNHXP])/;
		if ($1 >= $cmp)    {
                        $len=scalar(@CIGAR);
			$extract=$1 * -1;
			$base=substr(@_[9],$extract);
                        $qual=substr(@_[10],$extract);
			if (@_[1] & 0x10){
				$bitflag=165;
			} else	{
				$bitflag=133;
			}
                        print "@_[0]_$count\t$bitflag\t@_[2]\t@_[3]\t0\t*\t=\t@_[3]\t0\t$base\t$qual\n";
                        $togrep=$readlen-$1;
			$base=substr(@_[9],1,$togrep);
                        $qual=substr(@_[10],1,$togrep);
			$Ci=@_[5];
                        $Ci =~ s/$CIGAR[$#CIGAR]//;
                        print "@_[0]_$count\t73\t@_[2]\t@_[3]\t@_[4]\t$Ci\t@_[6]\t@_[3]\t0\t$base\t$qual\n";
                }
                $count++;
        }

}

	


