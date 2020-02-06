$head=qq{<HTML>
<HEAD>
<TITLE>Horizontal Gene Transfer Integration Site results</TITLE>
<META NAME="description" CONTENT="result">
<META NAME="keywords" CONTENT="result">
<META NAME="resource-type" CONTENT="document">
<META NAME="distribution" CONTENT="global">
<META HTTP-EQUIV="Content-Style-Type" CONTENT="text/css">
<!--[if IE]><script src="excanvas.js"></script><![endif]-->
<style type="text/css">
H1 { margin: 0 0 0 0; }
</style>
</HEAD>};

$appendix=qq{<H1><u><A NAME="Appendix"></A><BR>Appendix</H1></u>};
$candidate=qq{<BODY text="#000000" bgcolor="#FFFFFF" onload="">
<H1><BR>Candidate Integration Site</H1>
The following tables show Integration site candidates with gene annotation(<a href="#Integrationdetail">table description</a>).<BR>
<P><P><P><BR>};
$Integrationdetail=qq{<H1><A NAME="Integrationdetail"></A><BR>table description for Integration Site</H1><BR>
1. GeneName - HUGO gene name <BR>	
2. DistToGene - The reported distance for the nearest gene	<BR>
3. HumanChr - Integration Site Human chromosome 	<BR>
4. HIntegrationPos - exact Genomic location of Human chromosome for integration site	<BR> 
5. HChr	- Human chromosome for the region		<BR>
6. HStart - Genomic Start position of the cluster	<BR>
7. HPosition - Genomic End position of the cluster	<BR>
8. VirusID - Virus Name	<BR>
9. VStart - Genomic Start position of the viral cluster	<BR>
10. VirusEnd - Genomic End position of the viral cluster	<BR>
11. DiscordantReads - # of same read pair mapping to human cluster and Viral Cluster	<BR>
12. SoftclippingReads	- # of reads having Soft clip bases more than threshold<BR>
};

$primers=qq{
<H1><BR>Primer Design</H1>
The following tables show Primer design for the Integration Sites found(<a href="#Primerdetail">table description</a>).<BR>
<P><P><P><BR>};
$Primerdetail=qq{<H1><A NAME="Primerdetail"></A><BR>table description for Primer design</H1><BR>
1. Event - event ID<BR>
2. FWD_Sequence - left primer (before integration site within human genome for pair 1)/ (before integration site within viral genome for pair 2)<BR>
3. REV_Sequence - right primer (after integration site within viral genome for pair 1)/ (after integration site within human genome for pair 2)<BR>
4. FWD_TM - the predicated melting temperature of the left primer<BR>
5. REV_TM - the predicated melting temperature of the right primer<BR>
6. Product_Length - the length of the region being amplified<BR>
};
#### make a table for integration site detected
print $head;
print $candidate;
$tab=qq{<TABLE CELLPADDING=3 BORDER="1">};
$table=qq{</TABLE>
<P><P><P><BR>};
print $tab;
open FH, "$ARGV[0]" or die "";
##http://www.genecards.org/cgi-bin/carddisp.pl?gene=CCAT1&keywords=CCAT1
while(my $l = <FH>)	{
	chomp $l;
	@a=split(/\t/,$l);
	print "<TR>\n";	
	for($i=0;$i<=$#a;$i++)	{
		if ($. >1 && $i==0 && $a[$i] ne "NA")	{
		$entry=qq{<TD ALIGN="LEFT"><a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=$a[$i]&keywords=$a[$i]" target="_blank">$a[$i]</TD>};
		}else{
			$entry=qq{<TD ALIGN="LEFT">$a[$i]</TD>};
		}
		print "$entry\n";
	}
	print "</TR>\n";
}
print $table;
close FH;

print $primers;
print $tab;
open FH, "$ARGV[1]" or die "";
while(my $l = <FH>)     {
        chomp $l;
        @a=split(/\t/,$l);
        print "<TR>\n";
        for($i=0;$i<=$#a;$i++)  {
                $entry=qq{<TD ALIGN="LEFT">$a[$i]</TD>};
                print "$entry\n";
        }
        print "</TR>\n";
}
print $table;
close FH;

print qq{<link rel="stylesheet" href="fancybox/source/jquery.fancybox.css?v=2.0.6" type="text/css" media="screen" />
			<link rel="stylesheet" href="fancybox/source/helpers/jquery.fancybox-buttons.css?v=1.0.2" type="text/css" media="screen" />
			<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.7/jquery.min.js"></script>
			<script type="text/javascript" src="fancybox/source/jquery.fancybox.pack.js?v=2.0.6"></script>
			<script type="text/javascript" src="fancybox/source/helpers/jquery.fancybox-buttons.js?v=1.0.2"></script>

			<script type="text/javascript">
				\$(document).ready(function() {
					\$(".fancybox").fancybox({
								openEffect	: 'none',
								closeEffect	: 'none'
							});
				});
			</script>};


#### insert the Circos plot
### HUMAN
print "<H1>Circos plot showing Human Integration Sites</H1>\n";
$human_circos=qq{<a class="fancybox" rel="HGT-ID"  href="circosplot/human.circos.png" title="HGT-ID"> <img src="circosplot/human.circos.png" width="200"></a>};
print $human_circos;
### VIRAL
opendir DIR, "$ARGV[2]/circosplot" or die "cannot open dir $ARGV[2]/circosplot: $!";
my @file= readdir DIR;
print "<H1>Circos plot showing Viral Sites with overlayed coverage data</H1>\n";
#print "@file\n";
foreach $f (@file)	{
	if ($f ne 'human.circos.png' && $f ne '.' && $f ne '..')	{
		$name=$f;
		$name =~ s/.coverage.png//;
		
		$human_circos=qq{<a class="fancybox" rel="$name"  href="circosplot/$f" title="$name"> <img src="circosplot/$f" width="200"></a>};
		print $human_circos;
	}
}	
close DIR;
print $appendix;
print $Integrationdetail; 

print $Primerdetail;
print "</body>
</html>";