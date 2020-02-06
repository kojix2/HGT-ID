#!/usr/bin/perl

=head1 NAME
   ReadClustersdf.pl

=head1 SYNOPSIS

    USAGE: ReadCluster.pl --input input_sam_file --output output_prefix [--threshold 10000 --minClusterSize 4]

=head1 OPTIONS

B<--input,-i>
   Input file

B<--output,-o>
   output prefix

B<--window, -w>
    Window size

B<--minClusterSize, -m>
	Min size of cluster

B<--help,-h>
   This help message

=head1  DESCRIPTION


=head1  INPUT


=head1  OUTPUT


=head1  CONTACT
  baheti.saurabh@mayo.edu


==head1 EXAMPLE
   ReadCluster.pl --input=filename.sam --window=10000 --output=PREFIX

=cut

use strict;
#use warnings;
use Data::Dumper;
use DBI;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output|o=s',
                          'window|w=s',
						  'minClusterSize|m=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

my $r1_start = 0;
my $r2_start = 0;
my $r1_end = $r1_start + $options{window};
my $r2_end = $r2_start + $options{window};
my $r1_chr = "";
my $r2_chr = "";

my @cluster = ();

open (FHD, "<", $options{input}) or die "Cound not open file $options{input}\n";
open (INTER, ">", $options{output} . ".inter.sam") or die "Cound not open file $options{output}.inter.sam\n";

while (<FHD>){
	chomp $_;

	#skip processing lines starting with @ just print to output file.
	if ($_ =~ /^@/){
	#	print INTER $_."\n";
		next;
	}
	
	check_sequence($_);
}

close(FHD);
close(INTER);

exit(0);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = ("input", "output");

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	unless($options{window}) { $options{window} = 500; }
	unless($options{minClusterSize}) { $options{minClusterSize} = 4; }
}

#############################################################################
sub check_sequence {
	my $line = shift;

	my @data = split(/\t/, $line);
	#print "original=" . join ("\t",@data[0..7]) ."\n";
	## check if mates are within the window.
	#my $ii=inWindow($data[7], 2);
	#my $ii1=inWindow($data[3], 1);
	#print "$ii\t$ii1\n";
	#print "COORD=$r1_start\t$r2_start\t$r1_end\t$r2_end\t$r1_chr\t$r2_chr\n";
#	print "Read\t$line\n";
#	my $len=scalar(@cluster);
#	print "size:$len\n";
	#print "$data[3]\t$data[7]\t$data[2]\t$data[6]\n";
	next if ($data[2] eq "*" || $data[6] eq "*");
	if ((inWindow($data[3], 1)) && (inWindow($data[7], 2)) && ($r1_chr =~ /$data[2]/) && ($r2_chr =~ /$data[6]/)) {
		## if minClusterSize is reached output
		if (scalar(@cluster) >= $options{minClusterSize}) {
			#print "Print2File\t$line\n";
			print INTER $line."\n";
			#print "tothefile\n";

		} else {
			#print "ADD2cluster\t$line\n";
			push @cluster, $line;
			#print "add2cluster=" . join ("\t",@data[0..7]) ."\n";
		}
	} 
	else {
		#print "Else\t$line\n";
		#print "else=" . join ("\t",@data[0..7]) ."\n";
		if (scalar(@cluster) >= $options{minClusterSize}) {
			#print "dumping cluster\n";
			#print "@cluster\n";
			#my $len=scalar(@cluster);
			#print "size before dumping:$len\n";
			foreach (@cluster){
				print  INTER $_."\n";
			}
		}

		@cluster = ();
		$r1_start = $data[3] - $options{window};
		$r2_start = $data[7] - $options{window};
		#$r1_start = $data[3];
		#$r2_start = $data[7];
		$r1_end = $data[3] + $options{window};
		$r2_end = $data[7] + $options{window};
		$r1_chr = $data[2];
		$r2_chr = $data[6];
		#print "COORD-AfterElse=$r1_start\t$r2_start\t$r1_end\t$r2_end\t$r1_chr\t$r2_chr\n";
		#<STDIN>;
	}
	#<STDIN>;
}

#############################################################################
sub inWindow {
	my $coord = shift;
	my $read = shift;

	my $start = 0;
	my $end = 0;

	if ($read == 1) {
		$start = $r1_start- $options{window};
		$end = $r1_end+ $options{window};
	} else {
		$start = $r2_start- $options{window};
		$end = $r2_end+ $options{window};
	}
	#print "InWindow=$coord\t$start\t$end\n";
	if (($coord >= $start) && ($coord <= $end)){
		#print "#InWindow=$coord\t$start\t$end\t1\n";
		return 1;
	} else { #print "InWindow=$coord\t$start\t$end\t0\n";
	return 0; }
}

