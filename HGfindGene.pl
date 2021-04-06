#!/usr/bin/perl -w
use strict;
###################################
#get consistant probes either all upregulated significantly 
#cross a set of data or down-regulated significantly cross 
#a set of data
###################################
use CGI;
my $cgi= new CGI;
my $genes = $cgi->param('geneList');
my @stats=$cgi->param('stats');
my $PATH=$cgi->param('path');
my $ann_file=$cgi->param('long');
my $acr_file=$cgi->param('syn');
my $genelink=$cgi->param('link');
#my $factor_file=$cgi->param('factor');
#my $usage="$0 genes stat files\n";
#die $usage unless @ARGV >=2;
#my ($genes, @stats) = @ARGV;
#my $PATH="/data1/hairgrowth/final_stat/";
#my $genes="IL6";
#$stats[0]="FUE2.Biopsy.Gray.Thin_vs_NotThin";
#my $ann_file="/data1/chip_annotation/mouse/MmHs.long";
#my $acr_file="/data1/chip_annotation/mouse/MmHs.syn";
################################################
# 0. basic parameters
################################################

#my %factor;
#fill_hash(\%factor, $factor_file); 

my %ann;
fill_hash(\%ann, $ann_file);

my %acr;
fill_hash(\%acr, $acr_file);

my %PA;
foreach my $p(keys %acr){
    my $a=$acr{$p};
    if (! defined $PA{$a}){
	$PA{$a}=$p;
    }else{
	$PA{$a}.="\t".$p;
    }
}

#my %PB;
#foreach my $p(keys %factor){
#    my $a=$factor{$p};
#    if (! defined $PB{$a}){
#	$PB{$a}=$p;
#    }else{
#	$PB{$a}.="\t".$p;
#    }
#}




###########################################
# 1. clean the genelist get the probe set
###########################################
my @tmplist = split /\s*\n\s*|\t/, lc($genes);
my @probe;
my @miss;
my %search;
foreach my $gene (@tmplist){
    if ($gene =~ /^\s*$/){next;}
    if (defined $PA{$gene}){
	my $ps=$PA{$gene};
	my @tmp=split /\t/, $ps;
	foreach my $p (@tmp){
	    push(@probe, $p);
	    $search{$p}=$gene;
	}	
    }elsif (defined $acr{$gene}){
	my $ps=$PA{$acr{$gene}};
	my @tmp=split /\t/, $ps;
	foreach my $p (@tmp){
	    push(@probe, $p);
	   $search{$p}=$gene; 
	}
    }elsif (defined $ann{$gene}){
	push(@probe, $gene);
	$search{$gene}=$gene;
    #}elsif(defined $factor{$gene}){
       	
	#my $ps=$PA{$factor{$gene}};

	#my @tmp=split /\t/, $ps;
	#foreach my $p (@tmp){
	 #   push(@probe, $p);
	  #  $search{$p}=$gene;
	#}
 #   }elsif (defined $PB{$gene}){
#	my $ps=$PB{$gene};
#	my @tmp=split /\t/, $ps;
#	foreach my $p (@tmp){
#	    push(@probe, $p);
#	    $search{$p}=$gene;
#	}

    }else{
	push(@miss, $gene);
    }   	   
}



#######################################
# 2. load stat, and retrieve data
#      link = 'http://iip.pg.com:10080/iip-web/geneSummary.do?&annotationdata=20,56,13,74,44,11,55,50,51,49&acronym=' + gene;

########################################
print "Content-type: text/html\n\n";
print "<HTML><HEAD><TITLE>Gene Seach Result</TITLE>
<script language=\"JavaScript\">
    function gsr(gene){
	link='http://www.genecards.org/cgi-bin/carddisp.pl?gene='+gene;
      var fs = window.open( link, 'popup', 'width=900,height=700,resizable=yes,scrollbars=yes' );
      fs.focus();
      return (true);
    }

   function fue(gene){
      #link = 'http://mvic-biotech.na.pg.com/projects/Healthy_Hair/GSS0218_FUE/Statistical Results/PlotPerProbe/' + gene + '.PNG';
link = 'http://mvic-biotech.na.pg.com/projects/Hair_Growth/GSS0264_Little_Hair/Statistical%20Results/Plots/Graph/' + gene + '.PNG';

      var fs = window.open( link, 'popup', 'width=900,height=700,resizable=yes,scrollbars=yes' );
      fs.focus();
      return (true);
    }

  function cruella(gene){
      link = 'http://mvic-biotech.na.pg.com/projects/Healthy_Hair/GSS0255_Cruella/Statistical Results/PlotPerGene/' + gene + '.PNG';
      var fs = window.open( link, 'popup', 'width=900,height=700,resizable=yes,scrollbars=yes' );
      fs.focus();
      return (true);
    }


</script>\n</HEAD><BODY>\n";
my $header=my_color_header();
print $header;
print "<TABLE border=1>\n";
print "<H1>Gene Search Result</H1>\n";  
print "<TR><TH>&nbsp;</TH><TH>&nbsp;</TH>";   
foreach my $st (@stats){
    print "<TH colspan=2>$st</TH>";
}
print "<TH>&nbsp;</TH><TH>&nbsp;</TH></TR>\n";
get_all_stat(\@stats, \@probe, \%search, \%acr, \%ann);
print "</TABLE>\n";
print "</BODY></HTML>\n";

sub get_all_stat{
    my ($Hstats, $Hprobe, $Hsearch, $Hacr, $Hann)=@_;
    my %fc;
    my %pval;
    print "<TR><TH>SearchTerm</TH><TH>Affy_ID</TH>";
    foreach my $st (@{$Hstats}){
	my $stfile=$PATH.$st;
	if ($st =~ /siRNA_TGF/){
	   print "<TH>Average TGFb % Inhibition</TH><TH>Average % Viability</TH>"; 
	}else{
	    print "<TH>Fold</TH><TH>Pval</TH>";
	}
	fill_stat(\%{$fc{$st}}, \%{$pval{$st}}, $stfile);
    }
    
    print "<TH>Symbol</TH><TH>Description</TH></TR>\n";
    my %check;
    foreach my $gene (@{$Hprobe}){
	if (!defined $check{$gene}){
	    $check{$gene}=1; #only print once
	    my @all;
	    foreach my $st (@{$Hstats}){
		my @number=get_number($gene, \%{$fc{$st}}, \%{$pval{$st}});
		push (@all, @number);
	    }
	    my @number=get_number($gene, $Hacr, $Hann);
	    push (@all, @number);
	    my $mylink=$genelink.$gene.".PNG";
	    my $searchlink="https://www.google.com/?gws_rd=ssl#q=".$$Hsearch{$gene};
	    my $lline="<TR><TD><a href=\"".$searchlink."\">".$$Hsearch{$gene}."</a></TD><TD><a href=\"".$mylink."\">".$gene."</a></TD>";
	    my $pp=0;
	    for(my $i=0; $i<@all; $i++){
		my $n=$all[$i];
		
		if ($i== (@all-2)){
		    $n="<a href=\"".$searchlink."\">$n</a>";
		    $lline.="<TD>$n</TD>";
		   # print "<TD>$n</TD>";
		}elsif($n =~ /\#\#\#FC(\S+)/){
		    my $FC=$1;
		    if ($FC ne "NA"){$pp=1;}
		    if ($FC >4){ 
			$lline.="<td style=\"background-color:f06464\">$FC</td>"; 
		    }elsif($FC>2){
			$lline.="<td style=\"background-color:f09292\">$FC</td>";
		    }elsif($FC>1.2){
			$lline.="<td style=\"background-color:f0c1c1\">$FC</td>";
		    }elsif($FC>0.8){
			$lline.="<td style=\"background-color:f0f0f0\">$FC</td>";
		    }elsif($FC>0.5){
			$lline.="<td style=\"background-color:c1c1f0\">$FC</td>";
		    }elsif($FC>0.25){
			$lline.="<td style=\"background-color:9292f0\">$FC</td>";
		    }else{
			$lline.="<td style=\"background-color:6464f0\">$FC</td>";
		    }
		}elsif($n =~ /\#\#\#P(\S+)/){
		    my $P=$1;
		    if ($P ne "NA"){$pp=1;}
		    if ($P <0.0001){
			$lline.="<td style=\"background-color:ff0000\">$P</td>"; 
		    }elsif($P<0.001){
			$lline.="<td style=\"background-color:ff9700\">$P</td>";
		    }elsif($P<0.01){
			$lline.="<td style=\"background-color:ffca00\">$P</td>";
		    }elsif($P<0.1){
			$lline.="<td style=\"background-color:fffd00\">$P</td>";
		    }else{
			$lline.="<td style=\"background-color:c0c0c0\">$P</td>";
		    } 
		}else{
		    $lline.="<TD>$n</TD>";
		}
	    }
	    $lline.="</TR>";
	    if ($pp ==1 ){print $lline,"\n";}
	}
    }    
    return;
}



####################
# support methods
####################


sub my_color_header{
    my $header=
"
<table>
<tr><td valign=top>
<table border = 1>
<tr><th>FOLD_CHANGE_COLOR_KEY</th></tr>
<tr><td style=\"background-color:#6464f0\">Fold change &lt 0.25 </td></tr>
<tr><td style=\"background-color:#9292f0\"> 0.25 =&lt Fold change &lt 0.5</td></tr>
<tr><td style=\"background-color:#c1c1f0\"> 0.5 =&lt Fold change &lt 0.8</td></tr>

<tr><td style=\"background-color:#f0c1c1\"> 1.2 =&lt Fold change &lt 2</td></tr>
<tr><td style=\"background-color:#f09292\"> 2 =&lt Fold change &lt 4</td></tr>
<tr><td style=\"background-color:#f06464\"> 4 &lt Fold change</td></tr>

<tr><td style=\"background-color:f0f0f0\"> Other </td></tr>
</table>
</td>
<td valign=top>
<table border = 1>
<tr><th>pvalue_COLOR_KEY</th></tr>
<tr><td style=\"background-color:#FF0000\">p-value &lt 1.0E-4 </td></tr>
<tr><td style=\"background-color:#FF9700\"> 1.0E-4 =&lt p-value &lt 1.0E-3</td></tr>

<tr><td style=\"background-color:#FFCA00\"> 1.0E-3 =&lt p-value &lt 0.01</td></tr>
<tr><td style=\"background-color:#FFfd00\"> 0.01 =&lt p-value &lt 0.1</td></tr>
<tr><td style=\"background-color:#C0C0C0\"> other</td></tr>
</table>
</td>
<tr></table>
";
    return $header;

}
sub fill_hash{
    my ($h, $file)=@_;
    open (X, "$file")||die "could not open $file\n";
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;
	if (@tmp<2){next;}
	$$h{lc($tmp[0])}=lc($tmp[1]);
    }
    close X;
    return;
}
sub fill_factor{
    my ($h, $file)=@_;
    open (X, "$file")||die "could not open $file\n";
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, lc($line);
	if (@tmp<2){next};
	for(my $i=1;$i<@tmp;$i++){ 
	    $$h{$tmp[$i]}=$tmp[0];
	}
    }
    close X;
    return;
}

sub fill_stat{
    my ($h_fc,$h_p, $file)=@_;

    open (X, "$file")||die "could not open $file\n";
    my $title=<X>;
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;
	if (@tmp<3){next;}
	my $ff = sprintf("%.4f", $tmp[2]);
	$$h_fc{lc($tmp[0])}="###FC".$ff;
	$ff = sprintf("%.4f", $tmp[1]);
	$$h_p{lc($tmp[0])}="###P".$ff;

    }
    close X;
    return;
}
sub get_number{
    my ($name, @h)=@_;
    my @result;
    foreach my $hh (@h){
	if (defined $$hh{$name}){
	    push ( @result,$$hh{$name});
	}else{
	    push ( @result, "NA");
	}
    }
    return @result;
}
