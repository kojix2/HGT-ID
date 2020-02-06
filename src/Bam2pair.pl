#!/usr/bin/perl

use Cwd;
use File::Basename;
#Usage
sub usage(){
	print "Usage: perl Bam2Pair.pl -b <BAM> -o <outfile>\n
		-winsize [10000]\tThe distance between mate pairs to be considered the same\n
		-min [1]\t\tThe minimum number of reads required to support an SV event\n
		-prefix need a random prefix so files with the same name don't get created\n\n"
		;
}
$bedtools=`which intersectBed`;
$samtools=`which samtools`;

if(!defined($bedtools)){die "BEDtools must be installed\n";}
if(!defined($samtools)){die "Samtools must be installed\n";}
use Getopt::Long;
#Declare variables
GetOptions(
	'b=s' => \$BAM_FILE,		#path to bam
	'out=s' => \$outfile,		#path to output
	'winsize=i' => \$winsize,
    'prefix=s' => \$prefix,
	'min=i' => \$minSupport,
	'q=s' => \$qual,
	'v' => \$verbose,
	't' => \$keep_temp
	);
if(!defined($BAM_FILE)){die "Must specify a BAM file!\n".usage();}
if(!defined($outfile)){die "Must specify an out filename!\n".usage();}
if(!defined($qual)){$qual=20}


$Filter_BAM=$BAM_FILE;

@bam=split("/",$Filter_BAM);
$Filter_BAM=@bam[@bam-1];
$Filter_BAM=~s/.bam/$prefix.bam/;
$Filter_sam=$Filter_BAM;
$Filter_sam=~s/.bam/.sam/;




print "\nLooking for Discordant read pairs (and Unmated reads) without soft-clips\n";

### keep only the reads which are mapped to different chromsome
$command=join("","samtools view -f 1 -F 1800 -q $qual ",$BAM_FILE,"  |  awk -F\'\\t\' \'{if ((\$9 == 0 || \$1 ~ /^@/) && \$7 !~ /=/ ) print \$0}' > ",$Filter_sam);

print "$command\n" if ($verbose);
system($command);
$path = dirname(__FILE__);
$Filter_cluster=$Filter_sam;
$Filter_cluster=~s/.sam/.cluster/;
$command=join("",$path,"/ReadCluster.pl -i=$Filter_sam -o=$Filter_cluster -m=$minSupport ");
if($verbose){print "\n$command\n"};	

system($command);

$result_pe=join("",$Filter_cluster,".out");
$command=join("","cat ",$Filter_cluster,".inter.sam | sort -n -k3,3n -k4,4n |perl -ane 'if(\@F[6]!~/=/){print \"HGT\\t\" . join(\"\\t\",\@F[2],\@F[3],\@F[6],\@F[7],\"\\n\")}else{print \"HGT\\t\" , join(\"\\t\",\@F[2],\@F[3],\@F[2],\@F[7],\"\\n\")}' | sort -n -k2,12n -k3,12n >",$result_pe);
if($verbose){print $command."\n"};
system($command);
$command=join("","cat ",$result_pe," | ",$path,"/cluster.pair.pl ",$winsize," |awk '(\$6 >=",$minSupport,")' | awk 'BEGIN {OFS=\"\t\"} {if(\$0 !~ /^chr/){print \$3,\$4,\$1,\$2,\$6} else {print \$1,\$2,\$3,\$4,\$6}}'  | sort -n -k1,1n -k2,2n > ", $outfile);

##awk 'BEGIN {OFS="\t"} {if($0 !~ /^chr/){print $3,$4,$1,$2,$6} else {print $1,$2,$3,$4,$6}}'  | sort -n -k1,1n -k2,2n
if($verbose){print $command."\n"};
system($command);
$filt1=join("",$Filter_cluster,".inter.sam");
if (! $keep_temp ) {unlink($Filter_sam,$filt1,$result_pe);}
