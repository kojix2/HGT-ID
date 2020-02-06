#! /usr/bin/env bash

if [[ $# != 3 ]]
then
	echo "<input file from HGT-ID><output directory><configuration file>"
	exit 1;
fi

source $3
output=$2
### input is $1
if [[ `cat $1  | wc -l` -eq 1 ]]
then
	printf "Event\tFWD_Sequence\tREV_Sequence\tFWD_TM\tREV_TM\tProduct_Length\tEvent\tFWD_Sequence\tREV_Sequence\tFWD_TM\tREV_TM\tProduct_Length\n" > $output/output.primers
	exit 0;
fi	
#set -x
sed '1d' $1 | cut -f3,4 > $output/tmp_1-HumanInt
sed '1d' $1 | cut -f8,9 > $output/tmp_2-VirStart
sed '1d' $1 | cut -f8,10 > $output/tmp_3-VirEnd

cat $output/tmp_1-HumanInt | perl -lane '$start=$F[1]+50;$end=$F[1]+300;print"$F[0]:$start-$end";' > $output/tmp_1-HumanRight
cat $output/tmp_1-HumanInt | perl -lane '$start=$F[1]-300;$end=$F[1]-50;if($start<0){print"$F[0]:0-$end"}else{print"$F[0]:$start-$end"};' > $output/tmp_1-HumanLeft
cat $output/tmp_3-VirEnd | perl -lane '$start=$F[1]-300;$end=$F[1]-50;if($start<0){print"$F[0]:0-$end"}else{print"$F[0]:$start-$end"};' > $output/tmp_3-VirEnd_Left
cat $output/tmp_2-VirStart | perl -lane '$start=$F[1]+50;$end=$F[1]+300;print"$F[0]:$start-$end";' > $output/tmp_2-VirStart_Right

cat $output/tmp_1-HumanLeft | while read i; do $SAMTOOLS/samtools faidx $VIRUS_HUMAN_database $i >> $output/tmp_1-HumanLeft-Seq; done
cat $output/tmp_1-HumanRight | while read i; do $SAMTOOLS/samtools faidx $VIRUS_HUMAN_database $i >> $output/tmp_1-HumanRight-Seq; done
cat $output/tmp_2-VirStart_Right | while read i; do $SAMTOOLS/samtools faidx $VIRUS_HUMAN_database $i >> $output/tmp_2-VirStart_RightSeq; done
cat $output/tmp_3-VirEnd_Left | while read i; do $SAMTOOLS/samtools faidx $VIRUS_HUMAN_database $i >> $output/tmp_3-VirEnd_LeftSeq; done

left=$(cat $output/tmp_1-HumanLeft-Seq | tr '\n' ' ' | tr '>' '\n' | sed '1d' | cut -d' ' -f2- | sed s/' '//g); right=$(cat $output/tmp_2-VirStart_RightSeq | tr '\n' ' ' | tr '>' '\n' | sed '1d' | cut -d' ' -f2- | sed s/' '//g); echo $left | tr ' ' '\n' | perl -lane 'print"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN$F[0]";' > $output/tmp_4-FWD-1; echo $right | tr ' ' '\n' | perl -lane 'print"$F[0]NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";' > $output/tmp_4-FWD-2; cat $output/tmp_4-FWD-1 | rev | cut -b1-250 | rev > $output/tmp_4-FWD-3; cat $output/tmp_4-FWD-2 | cut -b1-250 > $output/tmp_4-FWD-4; paste $output/tmp_4-FWD-3 $output/tmp_4-FWD-4 | sed s/'\t'//g > $output/tmp_4-Front
 left=$(cat $output/tmp_3-VirEnd_LeftSeq | tr '\n' ' ' | tr '>' '\n' | sed '1d' | cut -d' ' -f2- | sed s/' '//g); right=$(cat $output/tmp_1-HumanRight-Seq | tr '\n' ' ' | tr '>' '\n' | sed '1d' | cut -d' ' -f2- | sed s/' '//g); echo $left | tr ' ' '\n' | perl -lane 'print"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN$F[0]";' > $output/tmp_4-REV-1; echo $right | tr ' ' '\n' | perl -lane 'print"$F[0]NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";' > $output/tmp_4-REV-2; cat $output/tmp_4-REV-1 | rev | cut -b1-250 | rev > $output/tmp_4-REV-3; cat $output/tmp_4-REV-2 | cut -b1-250 > $output/tmp_4-REV-4; paste $output/tmp_4-REV-3 $output/tmp_4-REV-4 | sed s/'\t'//g > $output/tmp_4-Back

cat -n $output/tmp_4-Front | perl -lane 'print"PRIMER_SEQUENCE_ID=LeftEvent_$F[0]\nSEQUENCE=$F[1]\nTARGET=249,2\nPRIMER_FIRST_BASE_INDEX=1\nPRIMER_TASK=pick_pcr_primers\nPRIMER_NUM_RETURN=1\nPRIMER_TM_SANTALUCIA=1\nPRIMER_SALT_CORRECTIONS=1\nPRIMER_DIVALENT_CONC=2\nPRIMER_DNTP_CONC=0.2\nPRIMER_DNA_CONC=400\nPRIMER_MIN_TM=57\nPRIMER_MAX_TM=59\nPRIMER_OPT_TM=58\nPRIMER_PRODUCT_SIZE_RANGE=100-600\nPRIMER_EXPLAIN_FLAG=1\n=";' > $output/tmp_5-Front.bio
cat -n $output/tmp_4-Back | perl -lane 'print"PRIMER_SEQUENCE_ID=LeftEvent_$F[0]\nSEQUENCE=$F[1]\nTARGET=249,2\nPRIMER_FIRST_BASE_INDEX=1\nPRIMER_TASK=pick_pcr_primers\nPRIMER_NUM_RETURN=1\nPRIMER_TM_SANTALUCIA=1\nPRIMER_SALT_CORRECTIONS=1\nPRIMER_DIVALENT_CONC=2\nPRIMER_DNTP_CONC=0.2\nPRIMER_DNA_CONC=400\nPRIMER_MIN_TM=57\nPRIMER_MAX_TM=59\nPRIMER_OPT_TM=58\nPRIMER_PRODUCT_SIZE_RANGE=100-600\nPRIMER_EXPLAIN_FLAG=1\n=";' > $output/tmp_5-Back.bio

cat $output/tmp_5-Front.bio | $PRIMER3 > $output/tmp_6-Front.details
cat $output/tmp_5-Back.bio | $PRIMER3 > $output/tmp_6-Back.details

cat $output/tmp_6-Front.details | egrep 'PRIMER_SEQUENCE_ID|PRIMER_LEFT_SEQUENCE|PRIMER_RIGHT_SEQUENCE|PRIMER_LEFT=|PRIMER_RIGHT=|LEFT_TM=|RIGHT_TM=' | sed s/'PRIMER_SEQUENCE_ID='/'!'/g | cut -d'=' -f2- | tr '\n' '\t' | tr '!' '\n' | sed '1d' | tr ',' '\t' | perl -lane '$math=(($F[6]+$F[5])-$F[3]);print"$F[0]\t$F[1]\t$F[2]\t$F[7]\t$F[8]\t$math";' > $output/tmp_7-Front.primers
cat $output/tmp_6-Back.details | egrep 'PRIMER_SEQUENCE_ID|PRIMER_LEFT_SEQUENCE|PRIMER_RIGHT_SEQUENCE|PRIMER_LEFT=|PRIMER_RIGHT=|LEFT_TM=|RIGHT_TM=' | sed s/'PRIMER_SEQUENCE_ID='/'!'/g | cut -d'=' -f2- | tr '\n' '\t' | tr '!' '\n' | sed '1d' | tr ',' '\t' | perl -lane '$math=(($F[6]+$F[5])-$F[3]);print"$F[0]\t$F[1]\t$F[2]\t$F[7]\t$F[8]\t$math";' > $output/tmp_7-Back.primers

printf "Event\tFWD_Sequence\tREV_Sequence\tFWD_TM\tREV_TM\tProduct_Length\tEvent\tFWD_Sequence\tREV_Sequence\tFWD_TM\tREV_TM\tProduct_Length\n" > $output/output.primers

sed s/'Left'/'Right'/g $output/tmp_7-Back.primers | paste $output/tmp_7-Front.primers - | sed s/'\t\t\t\t\0'/'NA\tNA\tNA\tNA\tNA'/g >> $output/output.primers

rm $output/tmp_*
