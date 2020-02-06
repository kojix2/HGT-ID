#!/usr/bin/perl -w
#!/usr/local/bin/perl5 -w

# ----- Installer Modifiable Variables -------------------------------------
# You may wish to modify the following variables to suit
# your installation.

$MAINTAINER
    ='<a href="mailto:your-email@your-host">your-email@your-host</a>';

# This is the select box to which you should
# add options for mispriming libraries.  Modify the
# %SEQ_LIBRARY variable in primer3_www_results as
# well.
$SELECT_SEQ_LIBRARY =
"<select name=PRIMER_MISPRIMING_LIBRARY>\n"
. "<option selected>NONE\n"

# Put "<option>'s"designating more repeat libraries here, e.g.
# . "<option> HUMAN\n"
# . "<option> RODENT\n"

. "</select><br>\n";

# Documentation for this screen.
$DOC_URL = "primer3_www_help.html";

# URL that will process the data from this screen.
$PROCESS_INPUT_URL = 'primer3_www_results.cgi';

# If you make any substantial modifications give this code a new
# version designation.
$CGI_VERSION = "(primer3_www.cgi v 0.1 beta 1a)";
# 1a has the objective function weight for product
# size reduced to .05.

# ----- End Installer Modifiable Variables ---------------------------------

$COPYRIGHT = q{ 
 Copyright (c) 1996,1997,1998
        Whitehead Institute for Biomedical Research. All rights reserved.
<p>
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
<ol>
<li>      Redistributions must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the  documentation
and/or other materials provided with the distribution.  Redistributions of
source code must also reproduce this information in the source code itself.

<li>      If the program is modified, redistributions must include a notice
(in the same places as above) indicating that the redistributed program is
not identical to the version distributed by Whitehead Institute.

<li>      All advertising materials mentioning features or use of this
software  must display the following acknowledgment:
<blockquote><i>
        This product includes software developed by the
        Whitehead Institute for Biomedical Research.
</i></blockquote>
<li>      The name of the Whitehead Institute may not be used to endorse or
promote products derived from this software without specific prior written
permission.
</ol>

We also request that use of this software be cited in publications as 
<blockquote><i>
Steve Rozen, Helen J. Skaletsky (1998)  Primer3. Code available at
<a href="http://www-genome.wi.mit.edu/genome_software/other/primer3.html">
http://www-genome.wi.mit.edu/genome_software/other/primer3.html.
</a>
</i></blockquote>
THIS SOFTWARE IS PROVIDED BY THE WHITEHEAD INSTITUTE ``AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE WHITEHEAD INSTITUTE BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.
};

BEGIN{
    print "Content-type: text/html\n\n";
    # Ensure that errors go to the web browser.
    open(STDERR, ">&STDOUT");
    $| = 1;
    print '';
}

use Carp;

use CGI;
# The CGI module is available from 
# http://www.genome.wi.mit.edu/ftp/distribution/software/WWW/

main();

sub main {

  $LIGHT_COLOR="#CCCCFF";
  $C_ROW = "<tr bgcolor=$LIGHT_COLOR>";
  $OLIGO_INPUT_SIZE=30;
  $SOURCE_SEQUENCE_WIDTH = 3 * ($OLIGO_INPUT_SIZE + 1) + 2;
  $SMALL_TEXT=9;
  $VSMALL_TEXT=4;

  $PR_DEFAULT_PRODUCT_MIN_SIZE = 100;
  $PR_DEFAULT_PRODUCT_OPT_SIZE = 200;
  $PR_DEFAULT_PRODUCT_MAX_SIZE = 1000;


  $PRIMER_SEQUENCE_ID                       = "";
  $SEQUENCE                                 = "";
  $INCLUDED_REGION                          = "";
  $TARGET                                   = "";
  $EXCLUDED_REGION                          = "";
  $PRIMER_SEQUENCE_QUALITY                  = "";
  $PRIMER_LEFT_INPUT                        = "";
  $PRIMER_RIGHT_INPUT                       = "";

  $PRIMER_MAX_MISPRIMING                    = "12.00";
  $PRIMER_PAIR_MAX_MISPRIMING               = "24.00";
  $PRIMER_GC_CLAMP                          = "0";
  $PRIMER_OPT_SIZE                          = "20";
  $PRIMER_MIN_SIZE                          = "18";
  $PRIMER_MAX_SIZE                          = "27";
  $PRIMER_OPT_TM                            = "60.0";
  $PRIMER_MIN_TM                            = "57.0";
  $PRIMER_MAX_TM                            = "63.0";
  $PRIMER_MAX_DIFF_TM                       = "100.0";
  $PRIMER_MIN_GC                            = "20.0";
  $PRIMER_OPT_GC_PERCENT                    = "";
  $PRIMER_MAX_GC                            = "80.0";
  $PRIMER_SALT_CONC                         = "50.0";
  $PRIMER_DNA_CONC                          = "50.0";
  $PRIMER_NUM_NS_ACCEPTED                   = "0";
  $PRIMER_SELF_ANY                          = "8.00";
  $PRIMER_SELF_END                          = "3.00";
  $PRIMER_MAX_POLY_X                        = "5";
  $PRIMER_NUM_RETURN                        = "5";
  $PRIMER_FIRST_BASE_INDEX                  = "1";
  $PRIMER_MIN_QUALITY                       = "0";
  $PRIMER_MIN_END_QUALITY                   = "0";
  $PRIMER_QUALITY_RANGE_MIN                 = "0";
  $PRIMER_QUALITY_RANGE_MAX                 = "100";
  $PRIMER_INSIDE_PENALTY                    = "";
  $PRIMER_OUTSIDE_PENALTY                   = "0";
  $PR_DEFAULT_MAX_END_STABILITY             = "9.0";
  $PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION    = "";
  $PRIMER_INTERNAL_OLIGO_INPUT              = "";
  $PRIMER_INTERNAL_OLIGO_OPT_SIZE           = "20";
  $PRIMER_INTERNAL_OLIGO_MIN_SIZE           = "18";
  $PRIMER_INTERNAL_OLIGO_MAX_SIZE           = "27";
  $PRIMER_INTERNAL_OLIGO_OPT_TM             = "60.0";
  $PRIMER_INTERNAL_OLIGO_MIN_TM             = "57.0";
  $PRIMER_INTERNAL_OLIGO_MAX_TM             = "63.0";
  $PRIMER_INTERNAL_OLIGO_MIN_GC             = "20.0";
  $PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT     = "";
  $PRIMER_INTERNAL_OLIGO_MAX_GC             = "80.0";
  $PRIMER_INTERNAL_OLIGO_SALT_CONC          = "50.0";
  $PRIMER_INTERNAL_OLIGO_DNA_CONC           = "50.0";
  $PRIMER_INTERNAL_OLIGO_SELF_ANY           = "12.00";
  $PRIMER_INTERNAL_OLIGO_MAX_POLY_X         = "5";
  $PRIMER_INTERNAL_OLIGO_SELF_END           = "12.00";
  $PRIMER_INTERNAL_OLIGO_MAX_MISHYB         = "12.00";
  $PRIMER_INTERNAL_OLIGO_MIN_QUALITY        = "0";
  $PRIMER_INTERNAL_OLIGO_NUM_NS             = "0";
  
  # $PRIMER_PRODUCT_MIN_TM                    = "-1000000.0";
  # $PRIMER_PRODUCT_OPT_TM                    = "0.0";
  # $PRIMER_PRODUCT_MAX_TM                    = "1000000.0";
  
  $PRIMER_WT_TM_GT                          = "1.0";
  $PRIMER_WT_TM_LT                          = "1.0";
  $PRIMER_WT_SIZE_LT                        = "1.0";
  $PRIMER_WT_SIZE_GT                        = "1.0";
  $PRIMER_WT_GC_PERCENT_LT                  = "0.0";
  $PRIMER_WT_GC_PERCENT_GT                  = "0.0";
  $PRIMER_WT_COMPL_ANY                      = "0.0";
  $PRIMER_WT_COMPL_END                      = "0.0";
  $PRIMER_WT_NUM_NS                         = "0.0";
  $PRIMER_WT_REP_SIM                        = "0.0";
  $PRIMER_WT_SEQ_QUAL                       = "0.0";
  $PRIMER_WT_END_QUAL                       = "0.0";
  $PRIMER_WT_POS_PENALTY                    = "0.0";
  $PRIMER_WT_END_STABILITY                  = "0.0";
  
  $PRIMER_IO_WT_TM_GT                       = "1.0";
  $PRIMER_IO_WT_TM_LT                       = "1.0";
  $PRIMER_IO_WT_SIZE_LT                     = "1.0";
  $PRIMER_IO_WT_SIZE_GT                     = "1.0";
  $PRIMER_IO_WT_GC_PERCENT_LT               = "0.0";
  $PRIMER_IO_WT_GC_PERCENT_GT               = "0.0";
  $PRIMER_IO_WT_COMPL_ANY                   = "0.0";
  $PRIMER_IO_WT_NUM_NS                      = "0.0";
  $PRIMER_IO_WT_REP_SIM                     = "0.0";
  $PRIMER_IO_WT_SEQ_QUAL                    = "0.0";
  
  $PRIMER_PAIR_WT_PR_PENALTY                = "1.0";
  $PRIMER_PAIR_WT_IO_PENALTY                = "0.0";
  
  $PRIMER_PAIR_WT_DIFF_TM                   = "0.0";
  $PRIMER_PAIR_WT_COMPL_ANY                 = "0.0";
  $PRIMER_PAIR_WT_COMPL_END                 = "0.0";
  $PRIMER_PAIR_WT_PRODUCT_TM_LT             = "0.0";
  $PRIMER_PAIR_WT_PRODUCT_TM_GT             = "0.0";
  
  $PRIMER_PAIR_WT_PRODUCT_SIZE_GT           = "0.05";
  $PRIMER_PAIR_WT_PRODUCT_SIZE_LT           = "0.05";
  
  $PRIMER_PAIR_WT_REP_SIM                   = "0.0";

  $query = new CGI;

  input_screen($query);
}

sub input_screen {
    my ($query) = @_;
    print $query->start_html("Primer3 Input $CGI_VERSION");

    my $input_buttons = qq{
        <input type="submit" name="Pick Primers" value="Pick Primers">
        <input type="reset" value="Reset Form"><br>
        };
    
    print qq{

<form action="$PROCESS_INPUT_URL" 
      method="post" 
      enctype="x-www-form-urlencoded">
<table border=0 width="100%">
  $C_ROW
    <td rowspan=2> <font size="+3"> Primer3 Test Pre-Release </font> <br>
        <a name="top"> pick primers from a DNA sequence</a>
    <td> <strong><a href="#disclaimer">disclaimer</a></strong> <br>
    <td rowspan=2> <a href="mailto:primer3\@genome.wi.mit.edu">
         bugs? suggestions?</a>
    <td rowspan=2> <a href="http://www.genome.wi.mit.edu/genome_software/other/primer3.html">
      source code</a>
  $C_ROW <td> <tt>   </tt><strong><a href=$DOC_URL#cautions>cautions</a></strong>
</table>

Paste <a name="PRIMER_SEQUENCE_INPUT"> source sequence</a>
below (5'->3', string of ACGTNacgtn --
other letters treated as N -- numbers and blanks ignored).
FASTA format ok.  Please N-out undesirable sequence
(vector, ALUs, LINEs, etc.)  or use a
<a name=PRIMER_MISPRIMING_LIBRARY_INPUT href="$DOC_URL#PRIMER_MISPRIMING_LIBRARY">
Mispriming Library (repeat library):</a>
};

    print $SELECT_SEQ_LIBRARY;

    print qq{

<textarea name="SEQUENCE" rows=6 cols=$SOURCE_SEQUENCE_WIDTH value="$SEQUENCE"></textarea>

<table border=0>
$C_ROW
 <td> <input type="checkbox" name="MUST_XLATE_PICK_LEFT" value=1 checked>
      Pick left primer or use left primer below.
 <td> <input type="checkbox" name="MUST_XLATE_PICK_HYB_PROBE">
      Pick hybridization probe (internal oligo) or use oligo below.
 <td> <input type="checkbox" name="MUST_XLATE_PICK_RIGHT" value=1 checked>
      Pick right primer or use right primer
      below (5'->3' on opposite strand).
<tr>
  <td> <input type="text" size=$OLIGO_INPUT_SIZE
        name=PRIMER_LEFT_INPUT value=$PRIMER_LEFT_INPUT>
  <td> <input type="text"  size=$OLIGO_INPUT_SIZE
        name=PRIMER_INTERNAL_OLIGO_INPUT value=$PRIMER_INTERNAL_OLIGO_INPUT>
  <td> <input type="text"  size=$OLIGO_INPUT_SIZE
        name=PRIMER_RIGHT_INPUT value=$PRIMER_RIGHT_INPUT>
</table>

$input_buttons

<table border=0 width="100%">

<tr> 
  <td><a name=PRIMER_SEQUENCE_ID_INPUT href="$DOC_URL#PRIMER_SEQUENCE_ID">Sequence Id:</a>
  <td><input type="text" name=PRIMER_SEQUENCE_ID value=$PRIMER_SEQUENCE_ID>
  <td> A string to identify your output.

<tr> 
 <td><a name=TARGET_INPUT href="$DOC_URL#TARGET">Targets:</a>
 <td><input type="text" name=TARGET value=$TARGET>
 <td>E.g. 50,2 requires primers to surround the 2 bases at positions 50 and 51.
     Or mark the <a href="#PRIMER_SEQUENCE_INPUT">source sequence</a>
     with [ and ]: e.g. ...ATCT[CCCC]TCAT.. means
     that primers must flank the central CCCC.


<tr> 
 <td><a name=EXCLUDED_REGION_INPUT href="$DOC_URL#EXCLUDED_REGION">Excluded Regions:</a>
 <td><input type="text" name=EXCLUDED_REGION value=$EXCLUDED_REGION>
 <td>E.g. 401,7 68,3 forbids selection of primers in the 7 bases starting at 401 and the 3 bases at 68.
     Or mark the <a href="#PRIMER_SEQUENCE_INPUT">source sequence</a>
     with &lt; and &gt;: e.g. ...ATCT&lt;CCCC&gt;TCAT.. forbids
     primers in the central CCCC.

</table>


<table border=0>
<tr>
 <td><a name=PRIMER_PRODUCT_SIZE_INPUT href="$DOC_URL#PRIMER_PRODUCT_SIZE">Product Size</a>
 Min: <input type="text" size=$SMALL_TEXT name=MUST_XLATE_PRODUCT_MIN_SIZE value=$PR_DEFAULT_PRODUCT_MIN_SIZE>
 Opt: <input type="text" size=$SMALL_TEXT name=PRIMER_PRODUCT_OPT_SIZE value=$PR_DEFAULT_PRODUCT_OPT_SIZE>
 Max: <input type="text" size=$SMALL_TEXT name=MUST_XLATE_PRODUCT_MAX_SIZE value=$PR_DEFAULT_PRODUCT_MAX_SIZE>
</table>

<table border=0>
<tr>
 <td><a name=PRIMER_NUM_RETURN_INPUT href="$DOC_URL#PRIMER_NUM_RETURN">Number To Return:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_NUM_RETURN value=$PRIMER_NUM_RETURN>
 <td><a name=PRIMER_MAX_END_STABILITY_INPUT href="$DOC_URL#PRIMER_MAX_END_STABILITY">Max 3\' Stability:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_MAX_END_STABILITY value=$PR_DEFAULT_MAX_END_STABILITY>

<tr>
 <td> <a name=PRIMER_MAX_MISPRIMING_INPUT href="$DOC_URL#PRIMER_MAX_MISPRIMING">Max Mispriming:</a>
 <td> <input type="text" size=$VSMALL_TEXT name=PRIMER_MAX_MISPRIMING value=$PRIMER_MAX_MISPRIMING>
 <td> <a name=PRIMER_PAIR_MAX_MISPRIMING_INPUT href="$DOC_URL#PRIMER_PAIR_MAX_MISPRIMING">Pair Max Mispriming:</a>
 <td> <input type="text" size=$VSMALL_TEXT name=PRIMER_PAIR_MAX_MISPRIMING value=$PRIMER_PAIR_MAX_MISPRIMING>
</table>

$input_buttons

<h3> General Primer Picking Conditions </h3>

<table border=0>
<tr>
 <td><a name=PRIMER_OPT_SIZE_INPUT href="$DOC_URL#PRIMER_SIZE">Primer Size</a>
 <td>Min: <input type="text" size=$VSMALL_TEXT name=PRIMER_MIN_SIZE value=$PRIMER_MIN_SIZE>
 <td>Opt: <input type="text" size=$VSMALL_TEXT name=PRIMER_OPT_SIZE value=$PRIMER_OPT_SIZE>
 <td>Max: <input type="text" size=$VSMALL_TEXT name=PRIMER_MAX_SIZE value=$PRIMER_MAX_SIZE>

<tr>
 <td><a name=PRIMER_OPT_TM_INPUT href="$DOC_URL#PRIMER_TM">Primer Tm</a>
 <td>Min: <input type="text" size=$VSMALL_TEXT name=PRIMER_MIN_TM value=$PRIMER_MIN_TM>
 <td>Opt: <input type="text" size=$VSMALL_TEXT name=PRIMER_OPT_TM value=$PRIMER_OPT_TM>
 <td>Max: <input type="text" size=$VSMALL_TEXT name=PRIMER_MAX_TM value=$PRIMER_MAX_TM>
 <td><a name=PRIMER_MAX_DIFF_TM_INPUT href="$DOC_URL#PRIMER_MAX_DIFF_TM">Max Tm Difference:</a>
     <input type="text" size=$VSMALL_TEXT name=PRIMER_MAX_DIFF_TM value=$PRIMER_MAX_DIFF_TM>

<tr>
 <td><a name=PRIMER_PRODUCT_TM_INPUT href="$DOC_URL#PRIMER_PRODUCT_TM">Product Tm</a>
 <td>Min: <input type="text" size=$VSMALL_TEXT name=PRIMER_PRODUCT_MIN_TM value=''>
 <td>Opt: <input type="text" size=$VSMALL_TEXT name=PRIMER_PRODUCT_OPT_TM value=''>
 <td>Max: <input type="text" size=$VSMALL_TEXT name=PRIMER_PRODUCT_MAX_TM value=''>

<tr>
 <td><a name=PRIMER_GC_PERCENT_INPUT href="$DOC_URL#PRIMER_GC_PERCENT">Primer GC%</a>
 <td>Min: <input type="text" size=$VSMALL_TEXT name=PRIMER_MIN_GC value=$PRIMER_MIN_GC>
 <td>Opt: <input type="text" size=$VSMALL_TEXT name=PRIMER_OPT_GC_PERCENT value=$PRIMER_OPT_GC_PERCENT>
 <td>Max: <input type="text" size=$VSMALL_TEXT name=PRIMER_MAX_GC value=$PRIMER_MAX_GC>

</table>

<table border=0>

<tr>
 <td><a name=PRIMER_SELF_ANY_INPUT href="$DOC_URL#PRIMER_SELF_ANY">Max Self Complementarity:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_SELF_ANY value=$PRIMER_SELF_ANY>
 <td><a name=PRIMER_SELF_END_INPUT href="$DOC_URL#PRIMER_SELF_END">Max 3\' Self Complementarity:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_SELF_END value=$PRIMER_SELF_END>

<tr>
 <td><a name=PRIMER_NUM_NS_ACCEPTED_INPUT href="$DOC_URL#PRIMER_NUM_NS_ACCEPTED">Max \#N\'s:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_NUM_NS_ACCEPTED value=$PRIMER_NUM_NS_ACCEPTED>

 <td><a name=PRIMER_MAX_POLY_X_INPUT href="$DOC_URL#PRIMER_MAX_POLY_X">Max Poly-X:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_MAX_POLY_X value=$PRIMER_MAX_POLY_X>

<tr>
  <td><a name=PRIMER_INSIDE_PENALTY_INPUT href="$DOC_URL#PRIMER_INSIDE_PENALTY">Inside Target Penalty:</a>
  <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INSIDE_PENALTY value=$PRIMER_INSIDE_PENALTY>
  <td><a name=PRIMER_OUTSIDE_PENALTY_INPUT href="$DOC_URL#PRIMER_OUTSIDE_PENALTY">Outside Target Penalty:</a>
  <td><input type="text" size=$VSMALL_TEXT name=PRIMER_OUTSIDE_PENALTY value=$PRIMER_OUTSIDE_PENALTY>
  <td><a name=PRIMER_INSIDE_PENALTY_INPUT href="$DOC_URL#PRIMER_INSIDE_PENALTY">
      Set Inside Target Penalty to allow primers inside a target.
      </a>

<tr>
 <td><a name=PRIMER_FIRST_BASE_INDEX_INPUT href="$DOC_URL#PRIMER_FIRST_BASE_INDEX">First Base Index:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_FIRST_BASE_INDEX value=$PRIMER_FIRST_BASE_INDEX>
 <td><a name=PRIMER_GC_CLAMP_INPUT href="$DOC_URL#PRIMER_GC_CLAMP">CG Clamp:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_GC_CLAMP value=$PRIMER_GC_CLAMP>

<tr>
  <td><a name=PRIMER_SALT_CONC_INPUT href="$DOC_URL#PRIMER_SALT_CONC">Salt Concentration:</a>
  <td><input type="text" size=$VSMALL_TEXT name=PRIMER_SALT_CONC value=$PRIMER_SALT_CONC>

  <td><a name=PRIMER_DNA_CONC_INPUT href="$DOC_URL#PRIMER_DNA_CONC">Annealing Oligo Concentration:</a>
  <td><input type="text" size=$VSMALL_TEXT name=PRIMER_DNA_CONC value=$PRIMER_DNA_CONC>

  <td><a name=PRIMER_DNA_CONC_INPUT href="$DOC_URL#PRIMER_DNA_CONC">
      (Not the concentration of oligos in the reaction mix but of those
       annealing to template.)</a>

</table>

<table border=0>
<tr>
 <td><input
  type="checkbox" name=PRIMER_LIBERAL_BASE value=1 checked>
  <a name=PRIMER_LIBERAL_BASE_INPUT 
     href="$DOC_URL#PRIMER_LIBERAL_BASE">Liberal Base</a>
				   
 <td><input
  type="checkbox" name=MUST_XLATE_PRINT_INPUT value=1>
  <a name=SHOW_DEBUGGING_INPUT 
     href="$DOC_URL#SHOW_DEBUGGING">Show Debuging Info</a>
</table>

$input_buttons

<h3> Other Per-Sequence Inputs </h3>

<table>
<tr><td><a name=INCLUDED_REGION_INPUT href="$DOC_URL#INCLUDED_REGION">Included Region:</a>
<td><input type="text" name=INCLUDED_REGION value=$INCLUDED_REGION>
<td> E.g. 20,400: only pick primers in the 400 base region starting at position 20.
    Or use { and } in the 
    <a href="#PRIMER_SEQUENCE_INPUT">source sequence</a>
    to mark the beginning and end of the included region: e.g.
    in ATC{TTC...TCT}AT the included region is TTC...TCT.

<tr><td><a name=PRIMER_START_CODON_POSITION_INPUT href="$DOC_URL#PRIMER_START_CODON_POSITION">Start Codon Position:</a>
<td><input type="text" name=PRIMER_START_CODON_POSITION value=''>

</table>
};

  print_seq_quality_input();

print qq{
<h3>Objective Function Penalty Weights for Primers</h3>

<table border=0>
};

obj_fn_weight_lt_gt('WT_TM', $PRIMER_WT_TM_LT, $PRIMER_WT_TM_GT, 'Tm');
obj_fn_weight_lt_gt('WT_SIZE', $PRIMER_WT_SIZE_LT, $PRIMER_WT_SIZE_GT, 'Size');
obj_fn_weight_lt_gt('WT_GC_PERCENT', $PRIMER_WT_GC_PERCENT_LT, $PRIMER_WT_GC_PERCENT_GT, 'GC%');

print "</table>\n\n<table border=0>";

obj_fn_weight('WT_COMPL_ANY', $PRIMER_WT_COMPL_ANY, 'Self Complementarity');
obj_fn_weight('WT_COMPL_END', $PRIMER_WT_COMPL_END, '3\' Self Complementarity');
obj_fn_weight('WT_NUM_NS', $PRIMER_WT_NUM_NS, '#N\'s');
obj_fn_weight('WT_REP_SIM', $PRIMER_WT_REP_SIM, 'Mispriming');
obj_fn_weight('WT_SEQ_QUAL', $PRIMER_WT_SEQ_QUAL, 'Sequence Quality');
obj_fn_weight('WT_END_QUAL',     $PRIMER_WT_END_QUAL,      'End Sequence Quality');
obj_fn_weight('WT_POS_PENALTY',  $PRIMER_WT_POS_PENALTY,   'Position Penalty');
obj_fn_weight('WT_END_STABILITY',$PRIMER_WT_END_STABILITY, 'End Stability');

print qq{
</table>
<h3>Objective Function Penalty Weights for Primer <it>Pairs</it></h3>
<table border=0>
};

obj_fn_weight_lt_gt('PAIR_WT_PRODUCT_SIZE', $PRIMER_PAIR_WT_PRODUCT_SIZE_LT,
    $PRIMER_PAIR_WT_PRODUCT_SIZE_GT, 'Product Size');
obj_fn_weight_lt_gt('PAIR_WT_PRODUCT_TM', $PRIMER_PAIR_WT_PRODUCT_TM_LT,
    $PRIMER_PAIR_WT_PRODUCT_TM_GT, 'Product Tm');

print "</table>\n\n<table border=0>\n";

obj_fn_weight('PAIR_WT_DIFF_TM',   $PRIMER_PAIR_WT_DIFF_TM,   'Tm Difference');
obj_fn_weight('PAIR_WT_COMPL_ANY', $PRIMER_PAIR_WT_COMPL_ANY, 'Any Complementarity');
obj_fn_weight('PAIR_WT_COMPL_END', $PRIMER_PAIR_WT_COMPL_END, '3\' Complementarity');
obj_fn_weight('PAIR_WT_REP_SIM',   $PRIMER_PAIR_WT_REP_SIM,   'Pair Mispriming');
obj_fn_weight('PAIR_WT_PR_PENALTY',$PRIMER_PAIR_WT_PR_PENALTY, 'Primer Penalty Weight');
obj_fn_weight('PAIR_WT_IO_PENALTY',$PRIMER_PAIR_WT_IO_PENALTY, 'Hyb Oligo Penalty Weight');

print "</table>\n\n";

print qq{
$input_buttons

<h3><a name="Internal+Oligo+Per-Sequence+Inputs">Hyb Oligo (Internal Oligo) Per-Sequence Inputs</a></h3>

<table border=0>
<tr><td><a name=internal_oligo_generic_INPUT href="primer3_www_doc.html#internal_oligo_generic">Hyb Oligo Excluded Region:</a>
<td><input type="text" name=PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION value=$PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION>

</table>

<h3><a name="Internal+Oligo+Global+Parameters">Hyb Oligo (Internal Oligo) General Conditions</a></h3>

<table border=0>
<tr>
 <td><a name=PRIMER_INTERNAL_OLIGO_SIZE_INPUT href="$DOC_URL#PRIMER_SIZE">Hyb Oligo Size:</a>
 <td>Min <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MIN_SIZE value=$PRIMER_INTERNAL_OLIGO_MIN_SIZE>
 <td>Opt <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_OPT_SIZE value=$PRIMER_INTERNAL_OLIGO_OPT_SIZE>
 <td>Max <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MAX_SIZE value=$PRIMER_INTERNAL_OLIGO_MAX_SIZE>

<tr>
 <td><a name=PRIMER_OPT_TM_INPUT href="$DOC_URL#PRIMER_TM">Hyb Oligo Tm:</a>
 <td>Min <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MIN_TM value=$PRIMER_INTERNAL_OLIGO_MIN_TM>
 <td>Opt <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_OPT_TM value=$PRIMER_INTERNAL_OLIGO_OPT_TM>
 <td>Max <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MAX_TM value=$PRIMER_INTERNAL_OLIGO_MAX_TM>

<tr>
 <td><a name=PRIMER_INTERNAL_OLIGO_GC_INPUT href="$DOC_URL#PRIMER_GC">Hyb Oligo GC%</a>
 <td>Min: <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MIN_GC value=$PRIMER_INTERNAL_OLIGO_MIN_GC>
 <td>Opt: <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT value=$PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT>
 <td>Max: <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MAX_GC value=$PRIMER_INTERNAL_OLIGO_MAX_GC>
</table>

<table border=0>

<tr><td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">
    Hyb Oligo Self Complementarity:</a>
    <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_SELF_ANY
               value=$PRIMER_INTERNAL_OLIGO_SELF_ANY>
    <td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">
        Hyb Oligo Max 3\' Self Complementarity:</a>
    <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_SELF_END
               value=$PRIMER_INTERNAL_OLIGO_SELF_END>

<tr>
 <td><a name=PRIMER_INTERNAL_OLIGO_NUM_NS_INPUT
      href="$DOC_URL#internal_oligo_generic">Max #Ns:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_NUM_NS
            value=$PRIMER_INTERNAL_OLIGO_NUM_NS>
 <td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">Hyb Oligo Max Poly-X:</a>
<td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MAX_POLY_X value=$PRIMER_INTERNAL_OLIGO_MAX_POLY_X>

<tr>
  <td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">Hyb Oligo Mishyb Library:</a>
  <td><select name=PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY>
      <option selected>NONE
      </select>

  <td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">Hyb Oligo Max Mishyb:</a>
  <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MAX_MISHYB
             value=$PRIMER_INTERNAL_OLIGO_MAX_MISHYB>

<tr><td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">
        Hyb Oligo Min Sequence Quality:</a>
    <td><input type="text" size=$VSMALL_TEXT
         name=PRIMER_INTERNAL_OLIGO_MIN_QUALITY value=$PRIMER_INTERNAL_OLIGO_MIN_QUALITY>

<tr><td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">Hyb Oligo Salt Concentration:</a>
    <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_SALT_CONC
               value=$PRIMER_INTERNAL_OLIGO_SALT_CONC>
    <td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">Hyb Oligo DNA Concentration:</a>
    <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_DNA_CONC
               value=$PRIMER_INTERNAL_OLIGO_DNA_CONC>
</table>

$input_buttons

<h3>Objective Function Penalty Weights for Hyb Oligos (Internal Oligos)</h3>

<table border=0>
};

obj_fn_weight_lt_gt('IO_WT_TM',   $PRIMER_IO_WT_TM_LT,   $PRIMER_IO_WT_TM_GT,   'Hyb Oligo Tm');
obj_fn_weight_lt_gt('IO_WT_SIZE', $PRIMER_IO_WT_SIZE_LT, $PRIMER_IO_WT_SIZE_GT, 'Hyb Oligo Size');
obj_fn_weight_lt_gt('IO_WT_GC_PERCENT', $PRIMER_IO_WT_GC_PERCENT_LT, $PRIMER_IO_WT_GC_PERCENT_GT, 'Hyb Oligo GC%');

print "</table>\n\n<table border=0>\n";


obj_fn_weight('IO_WT_COMPL_ANY', $PRIMER_IO_WT_COMPL_ANY, 'Hyb Oligo Self Complementarity');
obj_fn_weight('IO_WT_NUM_NS', $PRIMER_IO_WT_NUM_NS, 'Hyb Oligo #N\'s');
obj_fn_weight('IO_WT_REP_SIM', $PRIMER_IO_WT_REP_SIM, 'Hyb Oligo Mishybing');
obj_fn_weight('IO_WT_SEQ_QUAL', $PRIMER_IO_WT_SEQ_QUAL, 'Hyb Oligo Sequence Quality');

print qq{
</table>

$input_buttons

</form>

<h2><a name="disclaimer">Copyright Notice and Disclaimer</a></h2>
$COPYRIGHT
<h2>Acknowledgments</h2>
The development of Primer3 and the Primer3
web site was funded by 
<a href="http://www.hhmi.org/">Howard Hughes Medical Institute</a>
and by the 
<a href="http://www.nih.gov/">National Institutes of Health,</a>
<a href="http://www.nhgri.nih.gov/">
National Human Genome Research Institute.</a>
under grants R01-HG00257
(to David C. Page) and P50-HG00098 (to Eric S. Lander).
<p>
We gratefully acknowledge the support of
<a href="http://www.dec.com/">
Digital Equipment Corporation,</a>
which provided the Alphas which were used for much
of the development of Primer3, and of 
Centerline Software, Inc.,
whose TestCenter memory-error, -leak, and test-coverage checker
we use regularly to discover and correct otherwise latent
errors in Primer3.
<p>
<hr>
Original design of this primer-picking web site
by <em>Richard Resnick</em>, who
also is an author of this site\'s documentation.
Primer3\'s design is heavily based on earlier
implementations of similar programs:
Primer 0.5 (<em>Steve Lincoln, Mark Daly, and Eric
S. Lander</em>) and Primer v2 (<em>Richard Resnick</em>).
<em>Lincoln Stein</em> championed the use of the Boulder-IO
format and the idea of
making the Primer3 engine a software component.
<p>
Web software provided by
<em><a href="http://www.genome.wi.mit.edu/~steve">Steve Rozen</a></em>
<a href="mailto:steve\@genome.wi.mit.edu">steve\@genome.wi.mit.edu</a>
and <a href="http://www.genome.wi.mit.edu">
Whitehead Institute/MIT Center for Genome Research.</a>
<hr>
This page is maintained by $MAINTAINER.
<pre>$CGI_VERSION</pre>
};
    print $query->end_html;
    print "\n";
}

sub obj_fn_weight_lt_gt {
    my ($pretag, $def_lt, $def_gt, $row_head, $more) = @_;
    # <td><a name="${pretag}_INPUT" href="$DOC_URL#${pretag}">$row_head</a>
    print qq{
  <tr>
    <td><a name="${pretag}_INPUT" href="$DOC_URL#generic_penalty_weights">$row_head</a>
    <td>Lt:
    <td><input type=text size=$VSMALL_TEXT name="PRIMER_${pretag}_LT" value=$def_lt>
    <td>Gt:
    <td><input type=text size=$VSMALL_TEXT name="PRIMER_${pretag}_GT" value=$def_gt>
    };
    print "<td> $more\n" if $more;
}

sub obj_fn_weight {
    my ($tag, $def, $row_head, $more) = @_;
    #	<td><a name="${tag}_INPUT" href="$DOC_URL#$tag">$row_head</a>

    print qq{
	<tr>
	<td><a name="${tag}_INPUT" href="$DOC_URL#generic_penalty_weights">$row_head</a>
	<td><input type=text size=$VSMALL_TEXT name=PRIMER_$tag value=$def>
		};
    print "<td> $more\n" if $more;
}

sub print_seq_quality_input {
    print qq{
	<h3> <a name=PRIMER_SEQUENCE_QUALITY_INPUT href="$DOC_URL#PRIMER_SEQUENCE_QUALITY">Sequence Quality </a></h3>
	<textarea rows=2 cols=$SOURCE_SEQUENCE_WIDTH
	       name=PRIMER_SEQUENCE_QUALITY value=$PRIMER_SEQUENCE_QUALITY> </textarea>
        <table border=0>	    
<tr>
 <td><a name=PRIMER_MIN_QUALITY_INPUT href="$DOC_URL#PRIMER_MIN_QUALITY">Min Sequence Quality:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_MIN_QUALITY value=$PRIMER_MIN_QUALITY>
 <td><a name=PRIMER_MIN_END_QUALITY_INPUT href="$DOC_URL#PRIMER_MIN_END_QUALITY">Min End Sequence Quality:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_MIN_END_QUALITY value=$PRIMER_MIN_END_QUALITY>
 <td><a name=PRIMER_QUALITY_RANGE_MIN_INPUT href="$DOC_URL#PRIMER_QUALITY_RANGE_MIN">Sequence Quality Range Min:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_QUALITY_RANGE_MIN value=$PRIMER_QUALITY_RANGE_MIN>
 <td><a name=PRIMER_QUALITY_RANGE_MAX_INPUT href="$DOC_URL#PRIMER_QUALITY_RANGE_MAX">Sequence Quality Range Max:</a>
 <td><input type="text" size=$VSMALL_TEXT name=PRIMER_QUALITY_RANGE_MAX value=$PRIMER_QUALITY_RANGE_MAX>
       </table>
       };
}
