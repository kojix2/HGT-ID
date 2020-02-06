#use SeqComplexity;
use warnings;
$cCigar=12;


#ReadID  HChrom  HPos    HCiGar  Hseq    Virus   VPos    VCiGAR  VSeq
#HWI-ST590:243:C371BACXX:1:1101:1224:27485       4       180829216       51M     NGACCATCCTGGCCAACATGGTGAAACCCTGTCTCTACTAAAAATACAAAA
#     TorquetenomidivirusDNAisolate_AB303557  2631    32M19S  CCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAGAATGTTGAATATTNGCCC

while(<STDIN>){
	chomp;
	if($_=~/ReadID/){
		print $_,"\tVTemp\tVComplexity\n";
		next;
	}
	my @a=split(/\t/,$_);
	my $vcigar=$a[$cCigar];
	my $vseq=$a[$cCigar+1];
	my @cigarD=($vcigar=~m/(\d+)[SM]/g);
	#my $word=1; 
	my @cigarS=($vcigar=~m/[SM]/g);
	print $_,"\t";
	if(scalar(@cigarD) != scalar(@cigarS)){ 
		print "length of numbers differs length of characters\n";
		exit -1;
	}	
	#print "\nNEW: ",$vcigar,",",join(":",@cigarD),",",join(":",@cigarS),"\n";
	#print "VirusSeq: ",$vseq,"\n";
	my $max_cs=0;
	for (my $i=0;$i<scalar(@cigarD);$i++){
		#print "i=$i,cigarS=$cigarS[$i],cigarD=$cigarD[$i]\n";
		if($cigarS[$i] eq "S"){
			#print "S now\n";
			# clipped off
			$vseq=substr($vseq,$cigarD[$i]);
					
		}else{
		     if($cigarD[$i]>15){
			#print "M now\n";
			# match; calculate CompScore; clipped off;
			my $mseq=substr($vseq,0,$cigarD[$i]);
			my $rw=4;
			my $win=length($mseq);		
			#my @cs=();
			#for(my $n=1;$n<=4;$n++){
			#	push @cs, cl(\$mseq,\$win,\$n);
			#}
			#print "\t$mseq:",join(",",@cs);
			my @cs=cl(\$mseq,\$win,\$rw);
			$max_cs=($max_cs>$cs[0])?$max_cs:$cs[0];
			print " $mseq:",$cs[0];
			$vseq=substr($vseq,$cigarD[$i]);
		     }	
		}	
	}
	print "\t$max_cs\n";	
}





sub cl {
	my $k=4;
	my $seq    = shift;
	my $win    = shift;
	my $word   = shift;
	#print "inside sub ", $$seq," ",$$win," ",$$word,"\n";
	my $len    = length $$seq;
	my @values = ();
	#print STDERR "seq: $$seq\n";
	#print STDERR "length:", length($$seq),"\n";
	#print STDERR "win: $$win\n";
	#print STDERR "word: $$word\n";
	for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
		#print STDERR "p=$p\n";
		my $str     = substr ($$seq, $p, $$win);
		my $sum_vl  = 0;
		my $sum_vm  = 0;
		for (my $l = 1; $l <= $$word; $l++) {
			my $vl  = 0;
			my $vm  = 0;
			my $pot = pot($k, $l); 
			#print STDERR "k: $k; l: $l; pot: $pot\n";
			if ($pot < ($$win - $l + 1)) { $vm = $pot;           }
			else                         { $vm = $$win - $l + 1; }
			$sum_vm += $vm;
			#print STDERR "vm: $vm\n";
			my %elm = countWords($str, $l);
			foreach my $b (keys %elm) {
				next unless ($elm{$b} > 0);
				$vl++; 
			}
			$sum_vl += $vl;
		}
		my $r = 0;
		#print STDERR "sum_vl: $sum_vl\n";
		#print STDERR "sum_vm: $sum_vm\n";
		$r    = $sum_vl / $sum_vm if ($sum_vm > 0);
		push @values, $r;
	}
	#print "out values: @values\n";
	return @values;
}

=head3 pot

Function for calculate the exponential of a number.

Call: pot( $num, $exp ) NUMBER, NUMBER

Return: $res NUMBER

=cut

sub pot {
        my $num = shift @_;
        my $exp = shift @_;
        if ($num == 4) {
                if    ($exp == 1) { return       4; }
                elsif ($exp == 2) { return      16; }
                elsif ($exp == 3) { return      64; }
                elsif ($exp == 4) { return     256; }
                elsif ($exp == 5) { return    1024; }
                elsif ($exp == 6) { return    4096; }
                elsif ($exp == 7) { return   16384; }
                elsif ($exp == 8) { return   65536; }
                elsif ($exp == 9) { return  262144; }
                elsif ($exp ==10) { return 1048576; }
        }
        my $res = $num;
        for (my $i = 2; $i <= $exp; $i++) { $res *= $num; }
        return $res;
}

=head3 countWords

Function for count words in a sequence.

Call: countWords( $seq, $word ) STRING, NUMBER

Return: %results HASH (KEYS are the elements)

=cut

sub countWords {
        my $seq   = shift @_;
        my $word  = shift @_;
        my $len   = length $seq;
        my %count = ();
	my @alphabet = qw/A G C T/;
        # Init state when word == 1
        foreach my $b (@alphabet) {
                $count{$b} = 0;
        }

        for (my $i = 0; $i <= ($len - $word); $i++) {
                my $elem = substr ($seq, $i, $word);
                next if($elem =~ /[^ACGT]/);
                $count{$elem}++;
        }
        return %count;
}

