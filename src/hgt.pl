#!/usr/bin/perl

########################
### notes###
### 4/22
#used SM tag to read the sample name from the BAM header
##previously using ID tag
#### 9/13
### added the script to garner soft clip reads 
### it should be more than half of read length
### notes###
############################

#### Required bedtools & samtools to be in path
use Getopt::Long;
use strict;
use File::Basename;
use FindBin qw($Bin);
use Cwd;
use Cwd qw(chdir);
use File::Copy;
use File::Path qw(make_path remove_tree rmtree);
my ($BAMFILE,$SAMPLE,$OUTPUT,$CONFIG,$verbose,$debug);

#Declare variables
GetOptions(
	'b=s' => \$BAMFILE,
	's=s' => \$SAMPLE,
	'o=s' => \$OUTPUT,
	'c=s' => \$CONFIG,
	'v' => \$verbose, 
	'd' => \$debug,
	"help|h|?"	=> \&usage);

print "Usage = hgt.pl -b $BAMFILE -c $CONFIG \n\n";
sub usage {
	print "\nusage: hgt.pl [-sov] -b <BAMFILE> -c <CONFIG> \n";
	print "\t-s\t\tsample Name [Extracted from BAM header]\n";
	print "\t-o\t\tOutput Directory [cwd]\n";
	print "\t-v\t\tverbose flag [no]\n";
	print "\t-d\t\tdebug flag [no]\n";
	exit 1;
	}	
if(defined($BAMFILE)){$BAMFILE=$BAMFILE} else {print usage();die "Where is the BAM file?\n\n"}
if(defined($CONFIG)){$CONFIG=$CONFIG} else {print usage();die "Where is the configuration file?\n\n"}
##remove it
#rmtree("$OUTPUT");
if(defined($OUTPUT)){$OUTPUT=$OUTPUT;if (-d "$OUTPUT" || -e "$OUTPUT"){die "Output folder already exist\n";} elsif ( ! -d "$OUTPUT") {mkdir $OUTPUT}} else {$OUTPUT=getcwd;}
chdir ($OUTPUT) or die "cannot change: $!\n";
### creating the directory structure
my $logs=$OUTPUT ."/.logs";
my $autocode=$OUTPUT ."/.scripts"; 
mkdir $logs;
mkdir $autocode;
my $human_mapping=$OUTPUT ."/.human";
my $human_mapping_again=$OUTPUT ."/.human_again";
my $viral_mapping=$OUTPUT ."/.virus";
my $scoring=$OUTPUT ."/.scoring";
mkdir $human_mapping;
mkdir $human_mapping_again;
mkdir $viral_mapping;
mkdir $scoring;

#`bwa`;
### reading the configuration file
copy($CONFIG,"$OUTPUT/config.txt");
$CONFIG="$OUTPUT/config.txt";
my $config_vars=read_files_var($CONFIG);
my $SAMTOOLS=$config_vars->{'SAMTOOLS'};
my $BEDTOOLS=$config_vars->{'BEDTOOLS'};
my $CYTOBAND=$config_vars->{'HUMAN_CYTOBAND'};
my $BWA=$config_vars->{'BWA'};
my $PICARD=$config_vars->{'PICARD'};
my $FANCYBOX=$config_vars->{'FANCYBOX'};
my $RLIB=$config_vars->{'RLIB'};
### add tools to the path
$ENV{'PATH'} = $SAMTOOLS . ':' . $BEDTOOLS . ':' . $BWA . ':' . $PICARD . ':' . $ENV{'PATH'};
### get the path to all the scripts
#`bwa`;
my $SCRIPT_DIR=$Bin;
### references
my $USER_HUMAN_database=$config_vars->{'USER_HUMAN_database'};
my $HUMAN_database=$config_vars->{'HUMAN_database'};
my $HUMAN_database_Index=$config_vars->{'HUMAN_database_Index'};
my $VIRUS_database=$config_vars->{'VIRUS_database'};
my $VIRUS_database_Index=$config_vars->{'VIRUS_database_Index'};
my $VIRUS_HUMAN_database=$config_vars->{'VIRUS_HUMAN_database'};
my $VIRUS_HUMAN_database_Index=$config_vars->{'VIRUS_HUMAN_database_Index'};
my $REF_FLAT=$config_vars->{'REF_FLAT'};
my $MINRP=$config_vars->{'MINRP'};
my $MINSOFT=$config_vars->{'MINSOFT'};
### parameters 
my $THREADS=$config_vars->{'THREADS'};
my $SEQUENCE_COMPLEXITY=$config_vars->{'SEQUENCE_COMPLEXITY'};
my $DEPTH=$config_vars->{'DEPTH'};

### make sure the index file is available for the BAM file
my ($fn,$pathname) = fileparse($BAMFILE,".bam");
my $index=`ls $pathname/$fn*bai|head -1`;
if(!$index){die "\n\nERROR: you need to index your BAM file uisng samtools\n\n"}
### get current time
print "Start Time : " . &spGetCurDateTime() . "\n";
my $now = time;
my $command = "";
my $OUTNAME ="";
if ( ! $SAMPLE)	{
	$SAMPLE=`samtools view -H $BAMFILE |awk '{if(\$1~/^\@RG/){for(i;i<=NF;i++){if (\$i ~ /^SM:/){col=i}};sub("SM:","",\$col);print \$col}}' | head -1`;
	$SAMPLE=~s/\n//g;
}
my $libsize=`samtools view -q 20 -f2 $BAMFILE | awk '\$6 !~ /N/' | awk '\$7 ~ /=/' | cut -f9| awk '\$1<1000' | head -10000|awk '{if (\$1<0){\$1=-\$1}else{\$1=\$1} sum+=\$1;} END {printf(\"\%3.0f\\n",sum/NR)}'`;
chomp $libsize;
if ( $libsize =~ m/^\D/)        {
	$libsize=450;
}
my $readlength=`samtools view $BAMFILE | head -1 | awk '{print length(\$10)}'`;
chomp $readlength;
print "library size is inferred as : $libsize\n";
if($SAMPLE ne ""){$OUTNAME=$OUTPUT . "/" . $SAMPLE.".txt";}
else {$OUTNAME=$OUTNAME=$OUTPUT . "/output.txt";$SAMPLE="dummy";}
print "Writing results to $OUTNAME\n";
## to deal with TCGA sample/not mayo samples
my $chrflag=`samtools view -H $BAMFILE | grep '^\@SQ' | head -n1 | awk '{if(\$2 ~ /^SN:chr/) {print \"yes\"} else {print \"no\"}}'`;
chomp $chrflag;
##############################
### preprocessing done and calculations completed
################################################
print "Step1:\n";
print "extracting reads from human aligned bam for first pass mapping back to human...\n";
# $command=join ("","samtools view -u -b -f 8 -F 260 ",$BAMFILE," > ",$human_mapping,"/oneEndMapped.bam &");
### create a script to extract reads
open FH, ">$autocode/extract.h.sh" or die "can't open the script to write\n";
print FH "samtools view -u -b -f 8 -F 260 $BAMFILE > $human_mapping/oneEndMapped.bam \&\n";
print FH "pid=\$!\n";
print FH "samtools view -u -b -f 4 -F 264 $BAMFILE > $human_mapping/oneEndUnMapped.bam \&\n";
print FH "pid1=\$!\n";
### both reads are unmapped
print FH "samtools view -u -b -f 12 $BAMFILE > $scoring/UnMapped.bam \&\n";
print FH "pid4=\$!\n";


print FH "samtools view -F 8 -f 2 -F 4 -f 64 $BAMFILE  | awk '\$6 ~ /S/' | awk '\$6 !~ /I/ ' | awk '\$6 !~ /D/' | awk '\$6 ~ /S/' | perl -ane '\$len=length(\$F[9]);\$len=\$len/2;\@CIGAR = split(/([0-9]+[SMIDNHXP])/, \$F[5]);my \$hash;map { push(\@{\$hash->{\$2}}, \$1) if (/(\\d+)([SMIDNHXP])/) } \@CIGAR; foreach my \$softclip (\@{\$hash->{S}}) {if(\$softclip >= \$len){ print }}' > $human_mapping/read1 \&\n";
print FH "pid2=\$!\n";
print FH "samtools view -F 8 -f 2 -F 4 -f 128 $BAMFILE  | awk '\$6 ~ /S/' | awk '\$6 !~ /I/ ' | awk '\$6 !~ /D/' | awk '\$6 ~ /S/' | perl -ane '\$len=length(\$F[9]);\$len=\$len/2;\@CIGAR = split(/([0-9]+[SMIDNHXP])/, \$F[5]);my \$hash;map { push(\@{\$hash->{\$2}}, \$1) if (/(\\d+)([SMIDNHXP])/) } \@CIGAR; foreach my \$softclip (\@{\$hash->{S}}) {if(\$softclip >= \$len){ print }}'> $human_mapping/read2 \&\n";
print FH "pid3=\$!\n";
print FH "wait \$pid \$pid1 \$pid2 \$pid3 \$pid4\n";
print FH "##### get all the read IDs\n";
print FH "cat $human_mapping/read1 $human_mapping/read2  | cut -f1 | sort -T \$PWD | uniq > $human_mapping/IDS \n";
print FH "if [ -s $human_mapping/IDS ];\n then\n java -XX:ParallelGCThreads=$THREADS -Xmx6g -Xms3g -jar $PICARD/FilterSamReads.jar I=$BAMFILE WRITE_READS_FILES=false RLF=$human_mapping/IDS FILTER=includeReadList O=$human_mapping/softclip.bam VALIDATION_STRINGENCY=SILENT SORT_ORDER=queryname > $logs/0.FilterSamReads.log 2>&1;\n fi \n";
print FH "if [ -s $human_mapping/softclip.bam ];\n then\n samtools view $human_mapping/softclip.bam | perl $SCRIPT_DIR/splitReads.pl | samtools view -bt $HUMAN_database.fai - >  $human_mapping/soft.bam;\n fi \n";
close FH;
$command=join("","chmod 777 ",$autocode,"/extract.h.sh");
submit($command);
$command=join("","sh ",$autocode,"/extract.h.sh");
submit($command);

print "merging and sorting the Human BAM file...\n";
my $soft="$human_mapping/soft.bam";
if ( -e $soft)	{
$command=join ("","samtools merge -u -n ",$human_mapping,"/human.bam ",$human_mapping,"/oneEndMapped.bam ",$human_mapping,"/oneEndUnMapped.bam ",$human_mapping,"/soft.bam");	}
else{
	$command=join ("","samtools merge -u -n ",$human_mapping,"/human.bam ",$human_mapping,"/oneEndMapped.bam ",$human_mapping,"/oneEndUnMapped.bam ");
}
submit($command);
$command=join ("","samtools sort -\@ ",$THREADS," -n ",$human_mapping,"/human.bam ",$human_mapping,"/human.sort");
submit($command);
$command=join ("","samtools view -h ",$human_mapping,"/human.sort.bam | grep -v XF | samtools view -bS - > ",$human_mapping,"/tmp.bam");
submit($command);
$command=join("","mv ",$human_mapping,"/tmp.bam ",$human_mapping,"/human.sort.bam");
submit($command);
$command=join ("","perl ",$SCRIPT_DIR,"/keepreads.pl ",$human_mapping,"/human.sort.bam | sed '/^\\s*\$/d' | samtools view -u -bS - > ",$human_mapping,"/human.fix.sort.bam");
submit($command);
print "converting the BAM file to fastq...\n";
$command=join ("","java -XX:ParallelGCThreads=",$THREADS," -Xmx6g -Xms3g -jar ",$PICARD,"/SamToFastq.jar INPUT=",$human_mapping,"/human.fix.sort.bam FASTQ=",$human_mapping,"/read1.fq SECOND_END_FASTQ=",$human_mapping,"/read2.fq VALIDATION_STRINGENCY=SILENT > ",$logs,"/1.bam2fastq.log 2>&1");
submit($command);
#### get the FASTQ from BAM file for unmapped reads
$command=join ("","java -XX:ParallelGCThreads=",$THREADS," -Xmx6g -Xms3g -jar ",$PICARD,"/SamToFastq.jar INPUT=",$scoring,"/UnMapped.bam FASTQ=",$scoring,"/read1.fq SECOND_END_FASTQ=",$scoring,"/read2.fq VALIDATION_STRINGENCY=SILENT > ",$logs,"/1a.bam2fastq.log 2>&1");
submit($command);

#########################
#### completed step1
############################################################
#### map it again to human reference genome
print "Step2:\n";
print "mapping the reads back to human again...\n";
open FH, ">$autocode/run.h.bwa-mem.sh" or die "can't open the scipt to write\n";
print FH "bwa\n";
print FH "bwa mem -t $THREADS -M -R \'\@RG\\tID:$SAMPLE\\tSM:$SAMPLE\' $HUMAN_database_Index $human_mapping/read1.fq $human_mapping/read2.fq | samtools view -u -bS - >  $human_mapping_again/human.bam\n" ;
close FH;
$command=join ("","chmod -Rf 777 ",$autocode,"/run.h.bwa-mem.sh");
submit($command);
$command=join ("",$autocode,"/run.h.bwa-mem.sh > ",$logs,"/2.bwa-mem.log 2>&1");
submit($command);
print "extracting reads from again human aligned bam for viral mapping...\n";
open FH, ">$autocode/extract.h.again.sh" or die "can't open the script to write\n";
print FH "samtools view -u -b -f 8 -F 260 $human_mapping_again/human.bam > $human_mapping_again/oneEndMapped.bam \&\n";
print FH "pid=\$!\n";
print FH "samtools view -u -b -f 4 -F 264 $human_mapping_again/human.bam > $human_mapping_again/oneEndUnMapped.bam \&\n";
print FH "pid1=\$!\n";
print FH "wait \$pid \$pid1\n";
close FH;
$command=join("","chmod 777 ",$autocode,"/extract.h.again.sh");
submit($command);
$command=join("","sh ",$autocode,"/extract.h.again.sh");
submit($command);
print "merging and sorting the Human BAM file...\n";
$command=join ("","samtools merge -u -n ",$human_mapping_again,"/2_human.bam ",$human_mapping_again,"/oneEndMapped.bam ",$human_mapping_again,"/oneEndUnMapped.bam");
submit($command);
$command=join ("","samtools sort -@ ",$THREADS," -n ",$human_mapping_again,"/2_human.bam ",$human_mapping_again,"/human.sort");
submit($command);
$command=join ("","perl ",$SCRIPT_DIR,"/keepreads.pl ",$human_mapping_again,"/human.sort.bam | sed '/^\\s*\$/d' | samtools view -u -bS - > ",$human_mapping_again,"/human.fix.sort.bam");
submit($command);
print "converting the BAM file to fastq...\n";
$command=join ("","java -Xmx6g -Xms3g -jar ",$PICARD,"/SamToFastq.jar INPUT=",$human_mapping_again,"/human.fix.sort.bam FASTQ=",$human_mapping_again,"/read1.fq SECOND_END_FASTQ=",$human_mapping_again,"/read2.fq VALIDATION_STRINGENCY=SILENT > ",$logs,"/3.bam2fastq.log 2>&1");
submit($command);
### delete tmp files
if (! $debug){rmtree("$human_mapping");}
#########################
#### completed step2
############################################################
print "Step3:\n";
print  "mapping the reads to viral genome...\n";
open FH, ">$autocode/run.bwa-mem.sh" or die "can't open the scipt to write\n";
print FH "bwa\n";
print FH "bwa mem -t $THREADS -M -R \'\@RG\\tID:$SAMPLE\\tSM:$SAMPLE\' $VIRUS_database_Index $human_mapping_again/read1.fq $human_mapping_again/read2.fq | samtools view -u -bS - >  $viral_mapping/virus.bam \&\n";
print FH "pid=\$!\n";
print FH "bwa mem -t $THREADS -M -R \'\@RG\\tID:$SAMPLE\\tSM:$SAMPLE\' $VIRUS_database_Index $scoring/read1.fq $scoring/read2.fq | samtools view -u -bS - > $scoring/virus.fromUnmapped.bam \&\n";
print FH "pid1=\$!\n";
print FH "wait \$pid \$pid1\n";
close FH;
$command=join ("","chmod -Rf 777 ",$autocode,"/run.bwa-mem.sh");
submit($command);
$command=join ("",$autocode,"/run.bwa-mem.sh > ",$logs,"/4.bwa-mem.log 2>&1");
submit($command);
open FH, ">$autocode/extract.v.sh" or die "can't open the script to write\n";
print FH "samtools view -u -b -f 8 -F 260 $viral_mapping/virus.bam > $viral_mapping/oneEndMapped.bam \&\n";
print FH "pid=\$!\n";
print FH "samtools view -u -b -f 4 -F 264 $viral_mapping/virus.bam > $viral_mapping/oneEndUnMapped.bam \&\n";
print FH "pid1=\$!\n";
print FH "samtools sort -@ $THREADS $scoring/virus.fromUnmapped.bam $scoring/virus.fromUnmapped.sort \&\n";
print FH "pid2=\$!\n";
print FH "wait \$pid \$pid1 \$pid2\n";
print FH "samtools index $scoring/virus.fromUnmapped.sort.bam";
close FH;
$command=join("","chmod 777 ",$autocode,"/extract.v.sh");
submit($command);
$command=join("","sh ",$autocode,"/extract.v.sh");
submit($command);
$command=join ("","samtools merge -u -n ",$viral_mapping,"/2_virus.bam ",$viral_mapping,"/oneEndMapped.bam ",$viral_mapping,"/oneEndUnMapped.bam");
submit($command);
$command=join ("","samtools sort -@ ",$THREADS," -n ",$viral_mapping,"/2_virus.bam ",$viral_mapping,"/virus.sort");
submit($command);
$command=join ("","perl ",$SCRIPT_DIR,"/keepreads.pl ",$viral_mapping,"/virus.sort.bam | sed '/^\\s*\$/d' | samtools view -u -bS - > ",$viral_mapping,"/virus.fix.sort.bam");
submit($command);
#########################
#### completed step3
############################################################
print "Step4:\n";
print "finding the initial candidate list of HGT...\n";
#### create report for virus and human mapping for each HGT candidate 
$command=join ("","perl ",$SCRIPT_DIR,"/map.virus.human.pl ",$human_mapping_again,"/human.fix.sort.bam ",$viral_mapping,"/virus.fix.sort.bam ",$OUTPUT,"/HGT.candidates.txt ",$OUTPUT,"/VIRUS.HUMAN.sam > ",$logs,"/5.map.virus.human.log 2>&1");
submit($command);
$command=join ("","cat ",$OUTPUT,"/VIRUS.HUMAN.sam | sed '/^\$/d' | samtools view -bt ",$VIRUS_HUMAN_database,".fai - > ",$OUTPUT,"/VIRUS.HUMAN.bam");
submit($command);
$command=join ("","samtools sort -@ ",$THREADS," ",$OUTPUT,"/VIRUS.HUMAN.bam ",$OUTPUT,"/VIRUS.HUMAN.sort");
submit($command);
$command=join ("","samtools index ",$OUTPUT,"/VIRUS.HUMAN.sort.bam");
submit($command);
print "filtering candidate list uisng sequence complexity fiter ...\n";
$command=join ("","cat ",$OUTPUT,"/HGT.candidates.txt | perl ",$SCRIPT_DIR,"/filter.pl  | awk 'NR>1' | awk -v complex=$SEQUENCE_COMPLEXITY '\$NF>=complex' > ",$OUTPUT,"/HGT.candidates.flt.txt");
submit($command);
$command=join ("","java -Xmx6g -Xms3g -jar ",$PICARD,"/FilterSamReads.jar I=",$OUTPUT,"/VIRUS.HUMAN.sort.bam RLF=",$OUTPUT,"/HGT.candidates.flt.txt FILTER=includeReadList O=",$OUTPUT,"/VIRUS.HUMAN.flt.sort.bam WRITE_READS_FILES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=TRUE > ",$logs,"/6.FilterSamReads.log 2>&1");
submit($command);
#if ( ! $debug){unlink ("$OUTPUT/VIRUS.HUMAN.sort.reads");}
### get regions which are well covered
$command=join ("","samtools idxstats ",$OUTPUT,"/VIRUS.HUMAN.flt.sort.bam | awk -v depth=$DEPTH '\$NF+\$(NF-1)>depth' | cut -f1 > ",$OUTPUT,"/chromosomes2keep.txt");
submit($command); 
#########################
#### completed step4
############################################################
print "Step5:\n";
my $region2Capture=$readlength+$libsize;
print "creating regions to go back to original BAM file to extract reads for to find soft clipping...\n";
$command=join ("","bamToBed -i ",$OUTPUT,"/VIRUS.HUMAN.flt.sort.bam | sortBed -i stdin | mergeBed -i stdin | slopBed -i stdin -g ",$VIRUS_HUMAN_database,".fai -l ",$region2Capture," -r ",$region2Capture," | mergeBed -i stdin > ",$OUTPUT,"/regions.bed");
submit($command); 
my $chrs2keep=`cat $OUTPUT/chromosomes2keep.txt | tr \"\\n\" \"|\" | sed -e 's/\\(\.\*\\)./\\1/'`;
if ($chrflag eq "yes")	{
	$command=join("","cat ",$OUTPUT,"/regions.bed | grep -w -E '",$chrs2keep,"' > ",$OUTPUT,"/regions.2keep.bed");
}else{
	$command=join("","cat ",$OUTPUT,"/regions.bed | grep -w -E '",$chrs2keep,"' | sed -e 's/chrM/chrMT/g' | sed -e 's/chr//g' > ",$OUTPUT,"/regions.2keep.bed");
}
submit($command); 
if (! $debug) {unlink("$OUTPUT/HGT.candidates.txt","$OUTPUT/VIRUS.HUMAN.sam","$OUTPUT/VIRUS.HUMAN.bam","$OUTPUT/VIRUS.HUMAN.sort.bam","$OUTPUT/VIRUS.HUMAN.sort.bam.bai","$OUTPUT/HGT.candidates.flt.txt","$OUTPUT/VIRUS.HUMAN.sort.reads","$OUTPUT/chromosomes2keep.txt","$OUTPUT/VIRUS.HUMAN.flt.sort.reads","$OUTPUT/regions.bed");}

if (-z "$OUTPUT/regions.2keep.bed"){
	print "***********************************\n";
	print "Bad Luck !! NO HGT candidates found\n";
	print "***********************************\n";
	open OUT, ">$OUTNAME" or die "failed to open the file\n";
	print OUT "HumanChr\tHIntegrationPos\tHChr\tHStart\tHPosition\tVirusID\tVStart\tVirusEnd\tDiscordantReads\tSoftclippingReads\n";
	close OUT;
	if (! $debug) { rmtree("$viral_mapping","$human_mapping_again"); }
	exit 0;
}	
my $whichchr="";
my $chr2look=`cat $HUMAN_database.fai | cut -f1 |  tr \"\\n\" \"|\" | sed -e 's/\\(\.\*\\)./\\1/'`;
if ($chrflag == "yes")	{
	$whichchr=`cat $OUTPUT/regions.2keep.bed  | cut -f1 | sort | uniq | grep -w -E '$chr2look' | tr "\n"  " "| sed \'s/\\s*\$\/\/g\'`;
}else{
	$whichchr=`cat $OUTPUT/regions.2keep.bed  | awk '{print "chr"\$1}' | sort | uniq | grep -w -E '$chr2look'| tr "\n"  " "| sed \'s/\\s*\$\/\/g\'`;
}
chomp $whichchr;
print "go back to original BAM file to extract nearby reads...\n";

	
# $command=join ("","samtools view -b -f2 -F 256 -F 1024 -L ",$OUTPUT,"/regions.2keep.bed ",$BAMFILE," ",$whichchr," > ",$OUTPUT,"/human.org.bam");

$command=join ("","samtools view -b -F 256 -F 1024 -L ",$OUTPUT,"/regions.2keep.bed ",$BAMFILE," ",$whichchr," > ",$OUTPUT,"/human.org.bam");

submit($command); 
#exit;
#print "falg=$chrflag\n";
if ($chrflag eq "no" )	{
	$command=join ("","samtools view -h ",$OUTPUT,"/human.org.bam  | sed -e 's/SN:\\([0-9XYG]\\)/SN:chr\\1/' -e 's/SN:MT/SN:chrM/ -e 's/SN:chrM_rCRS/SN:chrM/' | samtools reheader - ",$OUTPUT,"/human.org.bam > ",$OUTPUT,"/human.org.tmp.bam");
	submit($command); 
	$command=join ("","mv ",$OUTPUT,"/human.org.tmp.bam ",$OUTPUT,"/human.org.bam");
	submit($command); 
}
$command=join ("","samtools view -b -f2 -F 256 -F 1024 -L ",$OUTPUT,"/regions.2keep.bed ",$human_mapping_again,"/human.bam > ",$OUTPUT,"/human.again.bam");
submit($command); 
$command=join ("","samtools view -b -f2 -F 256 -F 1024 -L ",$OUTPUT,"/regions.2keep.bed ",$viral_mapping,"/virus.bam > ",$OUTPUT,"/virus.org.bam");
submit($command); 
print "merge the BAM file to find the integration point...\n";
$command=join ("","samtools merge -f ",$OUTPUT,"/merged.bam ",$OUTPUT,"/VIRUS.HUMAN.flt.sort.bam ",$OUTPUT,"/human.org.bam ",$OUTPUT,"/human.again.bam ",$OUTPUT,"/virus.org.bam");
submit($command); 
$command=join ("","samtools sort -@ ",$THREADS," ",$OUTPUT,"/merged.bam ",$OUTPUT,"/",$SAMPLE,".forcalling");
submit($command); 
#$command=join ("","samtools index ",$OUTPUT,"/",$SAMPLE,".forcalling.bam");
#submit($command); 
#exit;
if (! $debug){rmtree("$human_mapping_again","$viral_mapping");}
if (! $debug){unlink ("$OUTPUT/merged.bam", "$OUTPUT/VIRUS.HUMAN.flt.sort.bam", "$OUTPUT/human.org.bam", "$OUTPUT/human.again.bam", "$OUTPUT/virus.org.bam","$OUTPUT/VIRUS.HUMAN.flt.sort.bai");}
#### filter the BAM file to remove the reads very close to 5KB centromere using the cytoband file

$command=join ("","zcat ",$CYTOBAND," | grep acen | awk '{if(\$4 ~ /^p/) {print \$1\"\\t\"\$2-50000\"\\t\"\$3} else {print \$1\"\\t\"\$2\"\\t\"\$3+50000}}' | intersectBed -abam ",$OUTPUT,"/",$SAMPLE,".forcalling.bam -b stdin -v > ",$OUTPUT,"/",$SAMPLE,".forcalling.cyto.bam" );
submit($command); 
$command=join ("","mv ",$OUTPUT,"/",$SAMPLE,".forcalling.cyto.bam ", $OUTPUT,"/",$SAMPLE,".forcalling.bam");
submit($command); 
$command=join ("","samtools index ",$OUTPUT,"/",$SAMPLE,".forcalling.bam");
submit($command); 

#########################
#### completed step5
############################################################
### find the integration points
print "Step6:\n";
print  "find the integration point...\n"; 

$command=join ("","perl ",$SCRIPT_DIR,"/integration.pl -b ",$OUTPUT,"/",$SAMPLE,".forcalling.bam -f ",$VIRUS_HUMAN_database, " -o ",$OUTNAME, " -d ",$libsize, " -m ",$MINRP," -l ",$MINSOFT," -r ",$HUMAN_database," -t -v ");
submit($command); 
#### filter the file to remove the reads very close to 5KB centromere using the cytoband file
$command=join("","cat ",$OUTNAME, " | awk '{if(NR==1){print \"#chr\tstart\tend\t\"\$0} else {print \$1\"\\t\"\$2\"\\t\"\$2\"\\t\"\$0}}' > ",$OUTNAME,".bed");
submit($command); 
$command=join ("","zcat ",$CYTOBAND," | grep acen | awk '{if(\$4 ~ /^p/) {print \$1\"\\t\"\$2-50000\"\\t\"\$3} else {print \$1\"\\t\"\$2\"\\t\"\$3+50000}}' | intersectBed -a ",$OUTNAME,".bed -b stdin -v -header | cut -f4- > ",$OUTNAME);
submit($command); 
if (! $debug){unlink ("$OUTNAME.bed");}
  
### map it to gene
print "create a gene bed file to annotate HGT...\n";
$command=join ("","zcat ",$REF_FLAT," | grep -v random | grep -v chrUn | grep -v hap | awk '{print \$3\"\\t\"\$5\"\\t\"\$6\"\\t\"\$1}' | sortBed -i stdin | perl ",$SCRIPT_DIR,"/uniqgene.pl > ",$OUTPUT,"/gene.bed");
submit($command); 
$command=join ("","cat ",$OUTNAME," | awk 'NR>1{print \$1\"\\t\"\$2\"\\t\"\$2}' | closestBed -t first -a stdin -d -b ",$OUTPUT,"/gene.bed  | awk '{print \$1\"\\t\"\$2\"\\t\"\$(NF-1)\"\\t\"\$NF}' > ",$OUTNAME,".gene.txt");
submit($command); 
$command=join ("","perl ",$SCRIPT_DIR,"/map.hgt.annot.pl ",$OUTNAME,".gene.txt ",$OUTNAME," | grep -v NA > ",$OUTNAME,".annot.txt");
submit($command); 
$command=join ("","mv ",$OUTNAME,".annot.txt ",$OUTNAME);
submit($command); 
if ( ! $debug){unlink("$OUTPUT/gene.bed","$OUTNAME.gene.txt","$OUTPUT/regions.2keep.bed");}
if ( ! $debug){rmtree("$human_mapping","$human_mapping_again","$viral_mapping");}

#### filter the results if it is both the human pairs

open VIRUSDB, "$VIRUS_database.fai" or die "can't open the VIRUS index file\n";
my %virus;
while(my $l=<VIRUSDB>)	{
	chomp $l;
	my @virus_contig=split(/\t/,$l);
	$virus{$virus_contig[0]}=1;
}
close VIRUSDB;
open OUT, "$OUTNAME" or die "can't open the output file\n";
open FLT, ">$OUTNAME.tmp.txt" or die "can't open the output folder\n";
while(my $l = <OUT>)	{
	print FLT $l if($.==1);
	chomp $l;
	my @call=split(/\t/,$l);
	if(defined $virus{$call[2]} || defined $virus{$call[7]})	{
		print FLT "$l\n";
	}	
}
close OUT;
close FLT;
$command=join ("","mv ",$OUTNAME,".tmp.txt ",$OUTNAME);
submit($command);


#### get the viral coverage 
my $coverage=$OUTPUT ."/.coverage";
mkdir $coverage;
my $count = `wc -l < $OUTNAME`;
sleep 10;

if ( $count > 1 )	{
	open FH, ">$autocode/coverage.sh" or die "can't open the script to write\n";
	print FH "for i in `cat $OUTNAME | awk 'NR>1' | cut -f8 | sort | uniq`\n";
	print FH "do\n";
	print FH "cat $VIRUS_database.fai | grep -w \$i | cut -f1,2 | perl $SCRIPT_DIR/coverage_split.pl | coverageBed -abam $viral_mapping/virus.fix.sort.bam -b stdin > $coverage/\$i.coverage.out \&\n";
	print FH "pid=\$!\n";
	print FH "cat $VIRUS_database.fai | grep -w \$i | cut -f1,2 | perl $SCRIPT_DIR/coverage_split.pl | coverageBed -abam $scoring/virus.fromUnmapped.sort.bam -b stdin > $coverage/\$i.proper.coverage.out \&\n";
	print FH "pid1=\$!\n";
	print FH "wait \$pid \$pid1\n";
	print FH "done\n";	
	close FH;
	$command=join("","chmod 777 ",$autocode,"/coverage.sh");
	submit($command);
	$command=join ("","sh ",$autocode,"/coverage.sh");
	submit($command);
	sleep 10;
	submit($command);
}
sleep 10;
### primer design
$command=join("","sh ",$SCRIPT_DIR,"/PrimerDesign_HGT_ID.sh ",$OUTNAME," ",$OUTPUT," ",$CONFIG);
submit($command);

#### create CIRCOS plots
my $plot=$OUTPUT ."/circosplot";
mkdir $plot;
$command=join("","Rscript ",$SCRIPT_DIR,"/circos.Rscript ",$OUTNAME," ",$OUTPUT," ",$plot," ",$RLIB);
submit($command);

#### get the scoring numbers
#get the Human Softclipped Reads
open FH, ">$autocode/scores.sh" or die "can't open the script to write\n";
print FH "echo \"HumanSoftClipping\" > $scoring/HSoft.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f5-7 | awk '{print \$1\":\"\$2-$readlength\"-\"\$3+$readlength}'`;\ndo\nsamtools view -f2 $OUTPUT/$SAMPLE.forcalling.bam \$i | awk '\$6 ~/S/' | awk '\$6 !~ /D/'  | awk '\$6 !~ /I/' | perl -ane 'my \@CIGAR = split(/([0-9]+[SMIDNHXP])/, \$F[5]); my \$hash; map { push(\@{\$hash->{\$2}}, \$1) if (/(\\d+)([SMIDNHXP])/) } \@CIGAR; \$count=0;foreach my \$softclip (\@{\$hash->{S}}) {if(\$softclip >= $MINSOFT){\$count++}}; print \"\$count\\n\"' | awk '\$1>0' | wc -l >> $scoring/HSoft.txt;\ndone\n";
print FH "echo \"ViralProperMappedReads\" > $scoring/VProper.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f8-10 | awk '{print \$1\":\"\$2\"-\"\$3}'`;\ndo\nsamtools view -f2 $scoring/virus.fromUnmapped.sort.bam \$i | cut -f1 | sort | uniq -c | wc -l >> $scoring/VProper.txt;\ndone\n";
print FH "echo \"ViralSoftClipping\" > $scoring/VSoft.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f8-10 | awk '{print \$1\":\"\$2\"-\"\$3}'`;\ndo\nsamtools view -f2 $scoring/virus.fromUnmapped.sort.bam \$i | awk '\$6 ~/S/' | awk '\$6 !~ /D/'  | awk '\$6 !~ /I/' | perl -ane 'my \@CIGAR = split(/([0-9]+[SMIDNHXP])/, \$F[5]); my \$hash; map { push(\@{\$hash->{\$2}}, \$1) if (/(\\d+)([SMIDNHXP])/) } \@CIGAR; \$count=0;foreach my \$softclip (\@{\$hash->{S}}) {if(\$softclip >= $MINSOFT){\$count++}}; print \"\$count\\n\"' | awk '\$1>0' | wc -l >> $scoring/VSoft.txt;\ndone\n";
print FH "echo \"HumanProperMappedReads\" > $scoring/HProper.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f5-7 | awk '{print \$1\":\"\$2-$readlength\"-\"\$3+$readlength}'`; do samtools view -f2 $OUTPUT/$SAMPLE.forcalling.bam \$i |cut -f1 | sort | uniq -c | wc -l  >> $scoring/HProper.txt;\ndone\n";
print FH "if [ $chrflag == \"yes\" ]\n";
print FH "then\n";
print FH "echo \"TotalCoverage\" > $scoring/TotalCoverageH.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f5-7 | awk '{print \$1\":\"\$2-$libsize\"-\"\$3+$libsize}'`; do samtools view $BAMFILE \$i | wc -l >> $scoring/TotalCoverageH.txt;\ndone\n";
print FH "else\n";
print FH "echo \"TotalCoverage\" > $scoring/TotalCoverageH.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f5-7 | awk '{print \$1\":\"\$2-$libsize\"-\"\$3+$libsize}' | sed -e 's/chr//g'`; do samtools view $BAMFILE \$i | wc -l >> $scoring/TotalCoverageH.txt;\ndone\n";
print FH "fi\n";
print FH "echo \"TotalCoverageViral\" > $scoring/TotalCoverageV.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f8-10 | awk '{print \$1\":\"\$2\"-\"\$3}'`; do samtools view $scoring/virus.fromUnmapped.sort.bam \$i | wc -l >> $scoring/TotalCoverageV.txt;\ndone\n";
print FH "paste $scoring/TotalCoverageV.txt $scoring/TotalCoverageH.txt | awk 'NR>1' | awk 'BEGIN{print \"Coverage\"} {print \$1+\$2}' > $scoring/TotalCoverage.txt\n";

print FH "cat $OUTNAME  | cut -f1-11 | paste - $scoring/HSoft.txt $scoring/HProper.txt $scoring/VProper.txt $scoring/VSoft.txt $scoring/TotalCoverage.txt  | awk '{if(NR==1){print \$0\"\\tScore\"} else {print \$0\"\\t\"(\$11+\$12+\$15-(((\$13+\$12)\*(\$14+\$15))/\$16))}}' > $OUTNAME.tmp.txt\n";
print FH "mv $OUTNAME.tmp.txt $OUTNAME\n";    
$command=join("","chmod 777 ",$autocode,"/scores.sh");
submit($command);
$command=join ("","sh ",$autocode,"/scores.sh");
submit($command);

if ( ! $debug){rmtree("$scoring");}

 
 
 
#### create an HTML page for the report and circos plot
if ( ! $debug){unlink("$viral_mapping/oneEndMapped.bam","$viral_mapping/2_virus.bam","$viral_mapping/oneEndUnMapped.bam","$viral_mapping/virus.bam","$viral_mapping/virus.fix.sort.bam","$viral_mapping/virus.sort.bam");} 
if ( ! $debug){rmtree("$human_mapping","$human_mapping_again","$viral_mapping");}

my $fancybox=$OUTPUT ."/fancybox";
mkdir $fancybox;
`cp -Rf $FANCYBOX/source $fancybox`;
$command=join("","perl ",$SCRIPT_DIR,"/create_html.pl ",$OUTNAME," ",$OUTPUT,"/output.primers ",$OUTPUT," > ",$OUTPUT,"/results.html");
submit($command);
   
print "Finish Time : " . &spGetCurDateTime() . "\n";
$now = time - $now;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60),
int($now % 60));

exit;
#########################
#### completed all the steps
############################################################

#######################
### subroutines
#######################
### parse the configuration file
sub read_files_var{
	my ($filename) = @_;
		open INFILE, "<$filename" or die "$filename $!\n";
	my %variables; 
	LINE:while (<INFILE>) {
	next LINE if /^#|^\s/;
	chomp;

	my ($var, $value) = split (/=/,$_);
	$value =~ tr/\"//d;
	$value =~ s/\s+$//;	
	$variables{$var}=$value;
	}
	return (\%variables);
}	
### get time
sub spGetCurDateTime {
	my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
	my $curDateTime = sprintf "%4d-%02d-%02d %02d:%02d:%02d",
	$year+1900, $mon+1, $mday, $hour, $min, $sec;
	return ($curDateTime);
}

sub submit{
	$command=shift;
	if($verbose){	print "$command\n";}
	system("$command");
}	
