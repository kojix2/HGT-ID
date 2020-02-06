use POSIX;
$windowsize=1000;
my %count;

while (<>){
    chomp;
    my ($chrstart, $start, $chrend, $end) = split /\t/;
    my $nstart = floor ($start/$windowsize);
    my $nend   = floor ($end/$windowsize);
    my $coord = "$start:$end";
	#print "$chrstart, $start, $chrend, $end, $coord\n";
    push @{$count{$chrstart}->{$nstart}->{$chrend}->{$nend}}, $coord;
}
print_groups (\%count);

sub print_groups {
    my ($rcount) = @_;
    my %count = %{$rcount};

    foreach my $chrstart (sort {$a<=>$b} keys %count) {
	foreach my $posstart (sort {$a<=>$b} keys %{$count{$chrstart}}) {
	    my %fcoord = %{$count{$chrstart}->{$posstart}};

	    foreach my $chrend (sort {$a<=>$b} keys %fcoord) {
		foreach my $posend (sort {$a<=>$b} keys %{$fcoord{$chrend}}){
		#  my $cpos = $fcoord{$chrend}->{$posend}->{start}[0];
		#	print "$chrstart,$posstart,$chrend,$posend\n";	
		#	print "$cpos\n";
			 my ($num, $avgx, $avgy) = calc_moments(@{$fcoord{$chrend}->{$posend}});
			 	
			print $chrstart."\t".$avgx."\t".$chrend."\t".$avgy ."\t".$num."\t" ;
			print "\n";
	    }
		}	
	}
	}
}	


sub calc_moments {
    my (@pos) = @_;
	#print "@pos\n";
    my ($num, $sumx, $sumy) = (0,0,0);
    foreach my $cpos (@pos) {
	@a=split(/:/,$cpos);
	$num++;
	$sumx+=$a[0];
	$sumy+=$a[1];
    }
    my $avgx = sprintf ("%d", $sumx/$num);
    my $avgy = sprintf ("%d", $sumy/$num);
	#print "$num\t$avgx,$avgy\n";
    return ($num, $avgx, $avgy);
}
