#!/bin/bash

#
#   Arguments check and usage information
#   Saurabh Baheti baheti.saurabh@mayo.edu
#   September 20 2016

usage()
{
cat << EOF
######
##      HGT-ID v1.0 installation script
##      Script Options:
##		-r	-		Reference genome FASTA file used to generate BAM file
##		-h	-		Display this usage/help text (No arg)
##		-v	-		verbose (No arg)
##
#################################################################################################################
##
## Authors:             Saurabh Baheti
## Creation Date:       September 20 2016
## Last Modified:       Januray 26 2017
##
## For questions, comments, or concerns, contact Saurabh Baheti(baheti.saurabh@mayo.edu)
##
#################################################################################################################
EOF
}

while getopts "r:vh" OPTION; do
  case $OPTION in
		v) verbose="YES" ;;
		h) usage
        exit ;;
        r) ref_genome=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG. See output file for usage." >&2
       usage
       exit ;;
    :) echo "Option -$OPTARG requires an argument. See output file for usage." >&2
       usage
       exit ;;
  esac
done

if [[ -z "$ref_genome"  && ! -s "$ref_genome" ]]
then
	echo "Must provide at least required options. See output file for usage." >&2
	usage
	exit 1;
fi

if [ "$verbose" ]; then set -x; fi

INSTALL_DIR=$(pwd)
BIN_DIR=$INSTALL_DIR/bin

echo -e "=========================================="
echo -e "HGT-ID PIPELINE (v1.0) installation script"
echo -e "=========================================="

SAMTOOLS="$BIN_DIR/samtools"
echo -e "Checking samtools directory...\n"
[ -d "$SAMTOOLS" ] || { echo -e "Error: directory '$SAMTOOLS' does not exist!" ; exit 1 ; }
echo -e "ok"

BEDTOOLS="$BIN_DIR/bedtools"
echo -e "Checking bedtools directory...\n"
[ -d "$BEDTOOLS" ] || { echo -e "Error: directory '$BEDTOOLS' does not exist!" ; exit 1 ; }
echo -e "ok"


BWA="$BIN_DIR/bwa"
echo -e "Checking bwa directory...\n"
[ -d "$BWA" ] || { echo -e "Error: directory '$BWA' does not exist!" ; exit 1 ; }
echo -e "ok"

PICARD="$BIN_DIR/picard"
echo -e "Checking picard directory...\n"
[ -d "$PICARD" ] || { echo -e "Error: directory '$PICARD' does not exist!" ; exit 1 ; }
echo -e "ok"

PRIMER3="$BIN_DIR/primer3"
echo -e "Checking primer3 directory...\n"
[ -d "$PRIMER3" ] || { echo -e "Error: directory '$PRIMER3' does not exist!" ; exit 1 ; }
echo -e "ok"

FANCYBOX="$BIN_DIR/fancybox"
echo -e "Checking fancybox directory...\n"
[ -d "$FANCYBOX" ] || { echo -e "Error: directory '$FANCYBOX' does not exist!" ; exit 1 ; }
echo -e "ok"

SCRIPTS="$INSTALL_DIR/src"
echo -e "Checking HGT-ID wrapper script directory...\n"
[ -d "$SCRIPTS" ] || { echo -e "Error: directory '$SCRIPTS' does not exist!" ; exit 1 ; }
echo -e "ok"

RESOURCES="$INSTALL_DIR/resources"
echo -e "Checking HGT-ID resources  directory...\n"
[ -d "$RESOURCES" ] || { echo -e "Error: directory '$RESOURCES' does not exist!" ; exit 1 ; }
echo -e "ok"

Rlibs="$INSTALL_DIR/Rlibs"
echo -e "Checking HGT-ID Rlibs  directory...\n"
[ -d "$Rlibs" ] || { echo -e "Error: directory '$Rlibs' does not exist!" ; exit 1 ; }
echo -e "ok"

#
# check for required executables
#

#### check for gcc complier
echo -e "Checking executables .....\n"
gcc=`gcc --version 2>&1`
if [[ ! `echo $gcc | grep gcc` ]]; then echo -e "no gcc compiler installed. Please install a gcc compiler "; exit 1;
else ver=`echo $gcc | awk '{print $3}'`;echo -e "gcc version: $ver                                      already installed, nothing to do ...\n";fi
echo -e "ok"

#### check and add Perl to env variable
perl=`which perl`
pp=`perl --version`
if [[ ! `echo $pp | grep Perl`  ]];then echo -e "No perl found, please install perl before using the tool ";exit 1;
else ver=`echo $pp | awk -F ',' '{print $2}'`;echo -e "perl version: $ver                                      already installed, nothing to do ...\n";fi
echo -e "ok"
perl=${perl%*/*}
#### java
java=`which java`
jj=`java -version 2>&1 | head -n 1 | awk -F '"' '{print $2}' | awk -F'.' 'BEGIN {OFS="."} {print $1,$2} '`
if [[ $(echo "$jj == 1.7" | bc) -ne 0 ]]; then echo -e "correct version of JAVA is installed";
else echo -e "1.7 version of JAVA is required to run the pipeline ....";exit 1; fi
if [[ ! $jj ]];then echo -e "No JAVA found, please install java before using the tool";exit 1;
else ver=` java -version 2>&1 | head -n 1| awk '{print $3}'`;echo -e "java version: $ver                                      already installed, nothing to do ...\n";fi
echo -e "ok"
java=${java%*/*}

#### Rexecutable
r_exe=`which R`
rr_version=`R --version 2>&1`
rr=`R --version 2>&1 |head -n 1 | awk '{print $3}' | awk -F'.' '{print $1}'`
if [[ $(echo "$rr == 3" | bc) -ne 0 ]]; then echo -e "correct version of R is installed";
else echo -e "version 3 or greater of R is required to run the pipeline ....";exit 1; fi
if [[ ! `echo $rr_version | grep 'R version'` ]];then echo -e "No R found, please install R before using the tool";exit 1;
else ver=`echo $rr_version | awk '{print $3}'`;echo -e "R version: $ver                                      already installed, nothing to do ...\n";fi
echo -e "ok"
r_exe=${r_exe%*/*}

### WGET execuatable 
wget=`which wget`
wg=`wget -V 2>&1`
if [[ ! `echo $wg | grep Wget` ]]
then
    echo "No commandline download tool found on your system. Please install wget on your machine and restart the installer"
	exit 1;
else
	ver=`echo $wg | awk '{print $3}'`
	echo "wget version: $ver                                      already installed, nothing to do ..."
fi


#
# Compile each package
#

echo -e "\nInstalling samtools...\n"
cd "$BIN_DIR/samtools"
`make --quiet > $BIN_DIR/samtools/log.txt 2>&1`
cd "$INSTALL_DIR"
echo -e "ok"

echo -e "Installing bedtools...\n"
cd "$BIN_DIR/bedtools"
bed=`make --quiet > $BIN_DIR/bedtools/log.txt 2>&1`
cd "$INSTALL_DIR"
echo -e "ok"
#rm log.txt

echo -e "Installing bwa...\n"
cd "$BIN_DIR/bwa"
`make --quiet > $BIN_DIR/bwa/log.txt 2>&1`
cd "$INSTALL_DIR"
echo -e "ok"

echo -e "Installing primer3...\n"
cd "$BIN_DIR/primer3/src"
`make --quiet > $BIN_DIR/primer3/log.txt 2>&1`
cd "$INSTALL_DIR"
echo -e "ok"

echo -e "Installing R packages ....\n"
cd $Rlibs
echo -e ".libPaths(\"$Rlibs\")" > $Rlibs/Install.R
echo -e "install.packages(\"circlize\", repos='http://cran.us.r-project.org')" >> $Rlibs/Install.R
$r_exe/Rscript $Rlibs/Install.R > $Rlibs/log.txt 2>&1
if [ $? -ne 0 ]
then
	echo -e "R packages failed to install ..."
	exit 1;
fi
	
cd ..

#clear

#### resources
echo -e "downloading and creating the indexes for the reference genome ....\n"
cp $ref_genome $RESOURCES/reference.fa
ref=$RESOURCES/reference.fa
echo -e "creating samtools index file ...\n"
$SAMTOOLS/samtools faidx $RESOURCES/reference.fa
echo -e "ok"

cd $RESOURCES
echo -e "downloading the Human reference genome .....\n"
# human reference genome
$wget -q ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta
cat Homo_sapiens_assembly19.fasta | awk '{if($0 ~ /^>/){print ">chr"$1} else {print}}' | sed -e 's/chr>/chr/g' | sed -e 's/chrMT/chrM/g' > human.fa
echo -e "indexing the human reference genome ....\n"
rm Homo_sapiens_assembly19.fasta
$SAMTOOLS/samtools faidx $RESOURCES/human.fa
echo -e "ok"
# viral reference genome
echo -e "downloading the viral reference genome .....\n"
$wget -q ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
zcat viral.1.1.genomic.fna.gz | awk -F'[|,]' '{if($0 ~/^>/){print ">"$5"_"$2} else {print}}' | tr " " "_" | sed -e 's/>_/>/g' > virus.fa
echo -e "indexing the viral reference genome ....\n"
rm viral.1.1.genomic.fna.gz
$SAMTOOLS/samtools faidx $RESOURCES/virus.fa
echo -e "ok"
### 
cat human.fa virus.fa > human_virus.fa
echo -e "indexing the human and viral reference genome ....\n"
$SAMTOOLS/samtools faidx $RESOURCES/human_virus.fa
echo -e "ok"
### BWA indexing the reference genome
echo -e "BWA indexing the reference genome ....\n"
#$BWA/bwa index -a bwtsw -p $RESOURCES/human_virus $RESOURCES/human_virus.fa & 
#pid=$!
$BWA/bwa index -a is -p $RESOURCES/virus $RESOURCES/virus.fa  & 
pid1=$!
$BWA/bwa index -a bwtsw -p $RESOURCES/human $RESOURCES/human.fa & 
pid2=$!
#wait $pid $pid1 $pid2
wait $pid1 $pid2
echo -e "ok"
# cytoband file
echo -e "downloading the cytoband file ...\n"
$wget -q http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz
# refFlat file
echo -e "downloading the refFlat file ...\n"
$wget -q http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refFlat.txt.gz


#blat=$HOME/bin/$MACHTYPE
echo -e "creating the config file to run your sample ....\n"
config=$INSTALL_DIR/src/config.txt
config=$( cat $config | tr "\n" "||")
config=${config//@INSTALL_DIR@/$INSTALL_DIR}
config=${config//@JAVA@/$java}
config=${config//@PERL@/$perl}
config=${config//@R_EXE@/$r_exe}
config=${config//@RESOURCES@/$RESOURCES}
echo $config | tr "||" "\n" > $INSTALL_DIR/config.txt
echo -e "ok"
cd ..

echo -e "HGT-ID succesfully is installed\n"
echo -e "The installation directory is $INSTALL_DIR\n"
echo -e "Add this path to your shell script environment file. (e.g, export HGT-ID=$INSTALL_DIR for bash) \n"






