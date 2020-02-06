#HWI-ST495:143193930:C2LRJACXX:1:2316:13994:79272        73      Humanherpesvirus7strainRK_AF037218      9622    0       37M63S  =       9622    0       CCTAACCCTAACCCTAACCCTAACCCAAACCCTAACCATAACACTGACCCTAACCGTAACATTGAATGTAATCATATGGAATAACCTAATGGTCTGGAAA       @<@DDDBDBA;BBFG9G@CF3<C3A3)11:1;B:))0*9BB###########################################################    NM:i:0  MD:Z:37 AS:i:37 XS:i:37 XA:Z:Humanherpesvirus7strainRK_AF037218,+152668,37M63S,0;  RG:Z:test
#HWI-ST495:143193930:C2LRJACXX:1:2316:13994:79272        133     Humanherpesvirus7strainRK_AF037218      9622    0       *       =       9622    0       CGAAAATTTCGTGTTATTCCATTGGTGGGGGTTTCGGTTAGGGTTACTTTTAAGAATAGTTCAGGCTGTTGCTTTTGGGATAGCTTTGCGTGTATCCTTT       ####################################################################################################    AS:i:0  XS:i:0  RG:Z:test-7B9A8B87


## read human bam 
## then read virus bam
#


open HUMAN , "samtools view $ARGV[0] | " or die "";
$read1=<HUMAN>;
$read2=<HUMAN>;
#chomp $read1;
#chomp $read2;
%calls;
print "parsing the HUMAN bam file ....\n";
while(defined $read1 && defined $read2)	{
	chomp $read1;
	chomp $read2;
	@R1=split(/\t/,$read1);
	@R2=split(/\t/,$read2);
	### find the orientation and read pair information
	if ($R1[1] & 0x10)	{
		$read1orient=-1;
	}else	{
		$read1orient=1;
	}
	####
	if ($R2[1] & 0x10)      {
			$read2orient=-1;
        }
	else   {
			$read2orient=1;
    }

	$value1=join("\t",@R1[0,2,3,5,9,10],1,$read1orient);
	#print "$value1\n";
	$value2=join("\t",@R2[0,2,3,5,9,10],2,$read2orient);	
	if ($R1[2] eq "\*" && $R2[2] eq "\*")	{ $read1=<HUMAN>; $read2=<HUMAN>;next;}
	if ($R1[5] eq "\*")	{
		$calls{$R1[0]}{1}=$value2;
		#print "$value2\n";
		#<STDIN>;
	}
	if ($R2[5] eq "\*")	{
		$calls{$R2[0]}{2}=$value1;
		#print "$value1\n";
		#<STDIN>;	
	}
	$read1=<HUMAN>;
	$read2=<HUMAN>;
	if ($. % 1000000 == 0)	{	print "$. lines for human read\n";} 

}
close HUMAN;
print "parsing the VIRUS bam file ....\n";
open VIRUS, "samtools view $ARGV[1]|" or die "";
$read1=<VIRUS>;
$read2=<VIRUS>;
chomp $read1;
chomp $read2;
open OUT , ">$ARGV[2]" or die "";
open SAM , ">$ARGV[3]" or die "";
print OUT "ReadID\tHChrom\tHPos\tHCiGar\tHseq\tHQual\tHR1/HR2\tHSeqReverse\tVR1/VR2\tVSeqReverse\tVirus\tVPos\tVCiGAR\tVSeq\tVQual\n";
while(defined $read1 && defined $read2) {
        @R1=split(/\t/,$read1);
        @R2=split(/\t/,$read2);
	 if ($R1[1] & 0x10)      {
                $read1orient=-1;
        }else   {
                $read1orient=1;
        }
        if ($R2[1] & 0x10)      {
                $read2orient=-1;
        }else   {
                $read2orient=1;
        }
	$value1=join("\t",1,$read1orient,@R1[2,3,5,9,10]);
	$value2=join("\t",2,$read2orient,@R2[2,3,5,9,10]);  

	# $value1=join("\t",@R1[2,3,5,9,10],1,$read1orient);
	# $value2=join("\t",@R2[2,3,5,9,10],2,$read2orient);
	

	if ($R1[5] eq "\*" )	{
	}
	else	{
		if (defined($calls{$R1[0]}{1}))	{
#			delete $calls{$R1[0]}{1};
			print OUT "$calls{$R1[0]}{1}\t$value1";
			#delete $calls{$R1[0]}{1};
			#### make a sam format
			($readid,$hchr,$hpos,$hcigar,$hseq,$hqual,$hwhichread,$hreverse)=split(/\t/,$calls{$R1[0]}{1});
			($vwhichread,$vreverse,$vchr,$vpos,$vcigar,$vseq,$vqual)=split(/\t/,$value1);
			if ($hwhichread == 1)	{
				($readid,$R1chr,$R1pos,$R1cigar,$R1seq,$R1qual,$R1whichread,$R1reverse)=split(/\t/,$calls{$R1[0]}{1});
				($R2whichread,$R2reverse,$R2chr,$R2pos,$R2cigar,$R2seq,$R2qual)=split(/\t/,$value1);	
			#	print "read1:$R1seq\n";
			}else	{
				($R1whichread,$R1reverse,$R1chr,$R1pos,$R1cigar,$R1seq,$R1qual)=split(/\t/,$value1);
				($readid,$R2chr,$R2pos,$R2cigar,$R2seq,$R2qual,$R2whichread,$R2reverse)=split(/\t/,$calls{$R1[0]}{1});

			}
			#print "read1: $R1seq\n";
			### read1 flag1
			### read2 flag2
			#print "$R1reverse \t $R2reverse\n";
			if ($R1reverse == -1 && $R2reverse == -1)	{
				$flag1=65+48;$flag2=129+48;
			}
			if ($R1reverse == -1 && $R2reverse != -1)	{
				$flag1=65+16;$flag2=129+32;
			}
			if ($R1reverse != -1 && $R2reverse == -1)       {
                $flag1=65+32;$flag2=129+16;
			}
			if ($R1reverse != -1 && $R2reverse != -1)       {
                $flag1=65;$flag2=129;
			}
			
			#print SAM join ("\t",$readid,$flag1,$R1chr,$R1pos,50,$R1cigar,$R2chr,$R2pos,0,$R1seq,$R1qual) . "\n";
			#print SAM join ("\t",$readid,$flag2,$R2chr,$R2pos,50,$R2cigar,$R1chr,$R1pos,0,$R2seq,$R2qual) . "\n";
			print SAM join ("\t",$readid,$flag1,$R1chr,$R1pos,50,$R1cigar,$R2chr,$R2pos,0,$R1seq,$R1qual) ."\n";
			print SAM join ("\t",$readid,$flag2,$R2chr,$R2pos,50,$R2cigar,$R1chr,$R1pos,0,$R2seq,$R2qual) . "\n";
			delete $calls{$R1[0]}{1};	



	}}
	if ($R2[5] eq "\*")      {
	}
	else	{
                if (defined($calls{$R2[0]}{2}))       {
#			delete $calls{$R1[0]}{2};
                        print OUT "$calls{$R2[0]}{2}\t$value2";
			    ($readid,$hchr,$hpos,$hcigar,$hseq,$hqual,$hwhichread,$hreverse)=split(/\t/,$calls{$R1[0]}{2});
                        ($vwhichread,$vreverse,$vchr,$vpos,$vcigar,$vseq,$vqual)=split(/\t/,$value2);
			if ($hwhichread == 1)   {
                                ($readid,$R1chr,$R1pos,$R1cigar,$R1seq,$R1qual,$R1whichread,$R1reverse)=split(/\t/,$calls{$R1[0]}{2});
                                ($R2whichread,$R2reverse,$R2chr,$R2pos,$R2cigar,$R2seq,$R2qual)=split(/\t/,$value2);
			 }else   {
                                ($R1whichread,$R1reverse,$R1chr,$R1pos,$R1cigar,$R1seq,$R1qual)=split(/\t/,$value2);
                                ($readid,$R2chr,$R2pos,$R2cigar,$R2seq,$R2qual,$R2whichread,$R2reverse)=split(/\t/,$calls{$R1[0]}{2});

                        }
			if ($R1reverse == -1 && $R2reverse == -1)	{
				$flag1=65+48;$flag2=129+48;
			}
			if ($R1reverse == -1 && $R2reverse != -1)	{
				$flag1=65+16;$flag2=129+32;
			}
			if ($R1reverse != -1 && $R2reverse == -1)       {
                $flag1=65+32;$flag2=129+16;
                        }
			
			#print SAM join ("\t",$readid,$flag1,$R1chr,$R1pos,50,$R1cigar,$R2chr,$R2pos,0,$R1seq,$R1qual) . "\n";
            #print SAM join ("\t",$readid,$flag2,$R2chr,$R2pos,50,$R2cigar,$R1chr,$R1pos,0,$R2seq,$R2qual) . "\n";

			print SAM join ("\t",$readid,$flag1,$R1chr,$R1pos,50,$R1cigar,$R2chr,$R2pos,0,$R1seq,$R1qual) ."\n" ;
            print SAM join ("\t",$readid,$flag2,$R2chr,$R2pos,50,$R2cigar,$R1chr,$R1pos,0,$R2seq,$R2qual)  ."\n" ;
	
			delete $calls{$R1[0]}{2};

        }	}
	$read1=<VIRUS>;
	$read2=<VIRUS>;
	if ($. % 10000 == 0)    {       print "$. lines for virus read\n";}
	if ( ! %calls)	{	last;}
	
}
close VIRUS;


close OUT;
				
	
	


