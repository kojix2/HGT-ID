#!/usr/bin/perl


#### Required bedtools & samtools to be in path
use Getopt::Long;
use strict;
use File::Basename;

my ($INPUT_BAM,$INPUT_FASTA,$OUTPUT_FILE,$minSoft,$dist_To_Soft,$bedtools,$samtools,$HUMAN_database);
my ($minRP,$MapQ,$minBQ);
my ($verbose,$keep_temp);

my $cmd = "";

#Declare variables
GetOptions(
	'b=s' => \$INPUT_BAM,
	'f=s' => \$INPUT_FASTA,
	'r=s' => \$HUMAN_database,
	'o:s' => \$OUTPUT_FILE,
	'm:i' => \$minRP,
	'l:i' => \$minSoft,
	't' => \$keep_temp,
	'd:i' => \$dist_To_Soft,
	'q:i' => \$MapQ,
	'v' => \$verbose, 
	"help|h|?"	=> \&usage);


if(defined($INPUT_BAM)){$INPUT_BAM=$INPUT_BAM} else {print usage();die "Where is the BAM file?\n\n"}
if(defined($INPUT_FASTA)){$INPUT_FASTA=$INPUT_FASTA} else {print usage();die "Where is the fasta file?\n\n"}
my ($fn,$pathname) = fileparse($INPUT_BAM,".bam");
my $index=`ls $pathname/$fn*bai|head -1`;
if(!$index){die "\n\nERROR: you need index your BAM file\n\n"}

### get current time
print "Start Time : " . &spGetCurDateTime() . "\n";
my $now = time;

if(defined($minSoft)){$minSoft=$minSoft} else {$minSoft=5}
if(defined($minRP)){$minRP=$minRP} else {$minRP=2}
if(defined($dist_To_Soft)){$dist_To_Soft=$dist_To_Soft} else {$dist_To_Soft=1000}
if(defined($MapQ)){$MapQ=$MapQ} else {$MapQ=20}

# adding and cheking for samtools and bedtools in the PATh
## check for bedtools and samtools in the path
$bedtools=`which intersectBed` ;
if(!defined($bedtools)){die "\nError:\n\tno bedtools. Please install bedtools and add to the path\n";}
#$samtools=`samtools 2>&1`;
$samtools=`which samtools`;
if($samtools !~ /(samtools)/i){die "\nError:\n\tno samtools. Please install samtools and add to the path\n";}


#Get sample name if available
###
my $SAMPLE_NAME="";
my $OUTNAME ="";
#$SAMPLE_NAME=`samtools view -f2 -H $INPUT_BAM|awk '{if(\$1~/^\@RG/){sub("ID:","",\$2);name=\$2;print name}}'|head -1`;
#$SAMPLE_NAME=~s/\n//g;
$SAMPLE_NAME=`samtools view -H $INPUT_BAM |awk '{if(\$1~/^\@RG/){for(i;i<=NF;i++){if (\$i ~ /^SM:/){col=i}};sub("SM:","",\$col);print \$col}}' | head -1`;
$SAMPLE_NAME=~s/\n//g;
if (!$OUTPUT_FILE){
	if($SAMPLE_NAME ne ""){$OUTNAME=$SAMPLE_NAME.".txt"}
	else {$OUTNAME="output.txt"}
}
else{$OUTNAME=$OUTPUT_FILE}

print "Writing results to $OUTNAME\n";

print "Usage = HGT.pl -l $minSoft -q $MapQ -d $dist_To_Soft -m $minRP -b $INPUT_BAM -f $INPUT_FASTA -o $OUTNAME \n\n";
sub usage {
	print "\nusage: HGT.pl [-cqlrmsd] -b <BAM> -f <Genome.fa> \n";
	print "\t-q\t\tMinimum mapping quality [20]\n";
	print "\t-l\t\tMinimum length of soft-clipped segment [5]\n";
	print "\t-m\t\tMinimum number of discordant read pairs [2]\n";
	print "\t-d\t\tMax distance between soft-clipped segments and discordant read pairs [1000]\n";
	print "\t-o\t\tOutput file name [output.txt]\n";
	print "\t-t\t\tPrint temp files for debugging [no|yes]\n";
	exit 1;
	}


#############################################################
# create temporary variable name
#############################################################
srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);
our $random_name=join "", map { ("a".."z")[rand 26] } 1..8;


my $tmp_name=join ("",$random_name,".tmp.bam");
my $random_file_sc = "";
my $command = "";

#############################################################
# Make sam file that has soft clipped reads
#############################################################
#give file a name
$random_file_sc=join ("",$random_name,".sc.sam");
$command=join ("","samtools view -q $MapQ -F 1024 $INPUT_BAM | awk '{OFS=\"\\t\"}{c=0;if(\$6~/S/){++c};if(c == 1){print}}' | perl -ane '\$TR=(\@F[10]=~tr/\#//);if(\$TR<2){print}' > ", $random_file_sc);
print "Making SAM file of soft-clipped reads\n";
if($verbose){	print "$command\n";}
system("$command");

#############################################################
# Find areas that have deep enough soft-clip coverage
print "Identifying soft-clipped regions that are at least $minSoft bp long \n";
open (FILE,"$random_file_sc")||die "Can't open soft-clipped sam file $random_file_sc\n";

my $tmpfile=join("",$random_file_sc,".sc.passfilter");
open (OUT,">$tmpfile")||die "Can't write files here!\n";

while(<FILE>){
	@_ = split(/\t/, $_);
	#### parse CIGAR string and create a hash of array of each operation
	my @CIGAR = split(/([0-9]+[SMIDNHXP])/, $_[5]);
	my $hash;
	map { push(@{$hash->{$2}}, $1) if (/(\d+)([SMIDNHXP])/) } @CIGAR;
	foreach my $softclip (@{$hash->{S}}) {
		if	($softclip > $minSoft){
			my $qual=$_[10];
			print OUT;
			last;
		}
	}
}
close FILE;
close OUT;

$command=join(" ","mv",$tmpfile,$random_file_sc);
if($verbose){	print "$command\n";}
system("$command");

#<STDIN>;
#########################################################
#Stack up SoftClips
#########################################################
my $random_file=join("",$random_name,".sc.direction.bed");
open (FILE,"$random_file_sc")|| die "Can't open sam file\n";
print "Calling sides of soft-clips\n";
open (TMPOUT,">$random_file")|| die "Can't create tmp file\n";

while (<FILE>){
	@_ = split(/\t/, $_);
	#### parse CIGAR string and create a hash of array of each operation
	my @CIGAR = split(/([0-9]+[SMIDNHXP])/, $_[5]);
	my $hash;
	map { push(@{$hash->{$2}}, $1) if (/(\d+)([SMIDNHXP])/) } @CIGAR;
	#### next if softclips on each end
	next if ($_[5] =~ /^[0-9]+S.*S$/);
	#### next softclip occurs in the middle
	next if ($_[5] =~ /^[0-9]+[^S][0-9].*S.+$/);
	my $softclip = $hash->{S}[0];
	my $end1 = 0;
	my $end2 = 0;
	my $softBases = "";
	my $right_corrected="";my $left_corrected="";

	if ($softclip > $minSoft) {
		####If the soft clip occurs at end of read and its on the minus strand, then it's a right clip
		if ($_[5] =~ /^.*S$/) {
			$end1=$_[3]+length($_[9])-$softclip-1;
			$end2=$end1+1;
			next if ($end1<0);
			#RIGHT clip on Minus
			$softBases=substr($_[9], length($_[9])-$softclip, length($_[9]));
			#print "$_[2],$end2,$softBases\n";
			$right_corrected=baseCheck($_[2],$end2,"right",$softBases);
            #print "$right_corrected\n";
			print TMPOUT "$right_corrected\n";
		} 
		else {
			#### Begins with S (left clip)
			$end1=$_[3]-$softclip;
			$end2=$end1+1;
			next if ($end1<0);
			$softBases=substr($_[9], 0,$softclip);
			$left_corrected=baseCheck($_[2],$end1,"left",$softBases);
			if(!$left_corrected){print "baseCheck($_[2],$end1,left,$softBases)\n";next}
			
            print TMPOUT "$left_corrected\n";
		}
	}
}
close FILE;
close TMPOUT;

## subroutine
sub baseCheck{
	my ($chrom,$pos,$direction,$softBases)=@_;
		#skip if position is less than 0, which is caused by MT DNA
		return if ($pos<0);
		my $exit="";

		while(!$exit){
		if($direction=~/right/){
			my $refBase=getSeq($chrom,$pos,$INPUT_FASTA);
			my $softBase=substr($softBases,0,1);
			if ($softBase !~ /$refBase/){
					my $value=join("\t",$chrom,$pos,$pos+1,join("|",$softBases,$direction));
					$exit="STOP";
					return $value;
			}
			else{
					$pos=$pos+1;
					$softBases=substr($softBases, 1,length($softBases));
			}
		 }
		else{
			my $refBase=getSeq($chrom,$pos+1,$INPUT_FASTA);
			my $softBase=substr($softBases,-1,1);
			if ($softBase !~ /$refBase/){
					$pos=$pos-1+length($softBases);
					my $value=join("\t",$chrom,$pos-1,$pos,join("|",$softBases,$direction));
					$exit="STOP";
					return $value;
			}
			else{
					$pos=$pos-1;
					$softBases=substr($softBases, 0, -1);
					#print "Trying again $softBases\n";
		   }
		}
	}
}
sub getSeq{
	my ($chr,$pos,$fasta)=@_;
	my @result=();
        if ($pos <0){print "$pos is not a valid position (likely caused by circular MT chromosome)\n";return;}
	@result = `samtools faidx $fasta $chr:$pos-$pos`;
	if($result[1]){chomp($result[1]);
	return uc($result[1]);
	}
	return("NA");
	#### after return will not be printed
	####print "RESULTS=@result\n";
}

#Remove SAM files to conserve space
if(! $keep_temp)	{unlink($random_file_sc);}

my $random_file_disc="$INPUT_BAM";
###
#
######################################################
# Transform Read pair groups into softclip equivalents
######################################################
my $v="";
my $t="";
print "Looking for discordant read pairs without requiring soft-clipping information\n";
use FindBin qw($Bin);
my $path=$Bin;
if($verbose){$v="-v"}
if($keep_temp){$t="-t"}
my $tmp_out=join("",$random_name,".out");
if ($keep_temp)	{
$command=join("","perl ",$path,"/Bam2pair.pl -b $random_file_disc -o $tmp_out -winsize $dist_To_Soft -min $minRP -prefix $random_name -q $MapQ $v $t");
}else	{
$command=join("","perl ",$path,"/Bam2pair.pl -b $random_file_disc -o $tmp_out -winsize $dist_To_Soft -min $minRP -prefix $random_name -q $MapQ $v");	
}
if($verbose){	print "$command\n"};
system("$command");
my $tmp_out_region=join("",$random_name,".regions.out");
$command=join("","perl ",$path,"/getregion.pl ",$tmp_out, " ",$dist_To_Soft," > ", $tmp_out_region);
if($verbose){	print "$command\n"};
system("$command");
#### remove the same genome calls
my $tmp_out_region1=join("",$random_name,".regions.1.out");
$command=join("","perl ",$path,"/filterSame.pl ",$tmp_out_region, " ",$HUMAN_database,".fai"," > ", $tmp_out_region1);
if($verbose){	print "$command\n"};
system("$command");
open (OUT,">$OUTNAME")|| die "Can't open file\n";
print OUT "HumanChr\tHIntegrationPos\tHChr\tHStart\tHPosition\tVirusID\tVStart\tVirusEnd\tDiscordantReads\tSoftclippingReads\n";

open (FILE,"$tmp_out_region1")|| die "Can't open sam file\n";
my $tmp_script=join("",$random_name,".sh");
while(my $l =<FILE>)	{
	next if ($. == 1);
	chomp $l;
	my @a=split(/\t/,$l);
	my $add=join(":",@a[0..2]);
	open (SCRIPT,">$tmp_script")|| die "Can't open script file\n";
	print SCRIPT "pos=`echo -e \"$a[0]\\t$a[1]\\t$a[2]\" | closestBed -a $random_file -b stdin -d -t first  | awk '\$NF< $dist_To_Soft &&  \$NF>0' | cut -f2 | sort | uniq -c | awk '{print \$2\"\\t\"\$1}' | sort -n -k2,2nr | head -1 | cut -f1`;\nreads=`echo -e \"$a[0]\\t$a[1]\\t$a[2]\" | closestBed -a $random_file -b stdin -d -t first  | awk '\$NF< $dist_To_Soft &&  \$NF>0' | wc -l`\nif [ \$pos ]; then echo -e \"$a[0]\\t\$pos\\t$l\\t\$reads\"; else pos=\$(echo \"scale=0; ($a[1]+$a[2])/2\" | bc -l); echo -e \"$a[0]\\t\$pos\\t$l\\t\$reads\";fi >>$OUTNAME";
	$command=join("","chmod 777 ",$tmp_script);
	if($verbose){	print "$command\n"};
	system($command);
	$command=join("","sh ",$tmp_script);
	if($verbose){	print "$command\n"};
	system($command);
	close SCRIPT;
	#<STDIN>;
}
close FILE;	

### filter the calls which are not HGT candidate

if(! $keep_temp)        {unlink($tmp_out_region,$tmp_out,$random_file,$tmp_script);}
print "Finish Time : " . &spGetCurDateTime() . "\n";
$now = time - $now;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60),
int($now % 60));

exit;

### sub routines
###############################################################################
#### to get time
sub spGetCurDateTime {
	my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
	my $curDateTime = sprintf "%4d-%02d-%02d %02d:%02d:%02d",
	$year+1900, $mon+1, $mday, $hour, $min, $sec;
	return ($curDateTime);
}
