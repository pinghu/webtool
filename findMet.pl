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
my $full_file=$cgi->param('long');

my $genelink=$cgi->param('link');
#my $factor_file=$cgi->param('factor');
#my $usage="$0 genes stat files\n";
#die $usage unless @ARGV >=2;
#my ($genes, @stats) = @ARGV;
#my $PATH="/home/ping/project/MicroSkin2014/2015/metabolon15/metabolon_stat/";
#my $genes="C00037";
#$stats[0]="M400.MlutS95M23XV2.PGAM715";
#my $full_file="/home/ping/project/MicroSkin2014/2015/metabolon15/metabolon14_15.xls";


###########################################
# 1. clean the genelist get the probe set
###########################################
my @tmplist = split /\s*\n\s*|\t/, lc($genes);
my @probe;
my @miss;

my %searchTerm;
foreach my $search (@tmplist){
    chomp $search;
    for ($search){
	s/\^s+//gi;
	s/\s+$//gi;
	s/\s+/ /gi;
    }
    $searchTerm{lc($search)}=0;
}

#my %searchID;
my %search;
my @ID;
my %ann;
open (ZZ, "$full_file")||die "could not open $full_file\n";
my $ttt=<ZZ>;
#chomp $ttt;
#for($ttt){s/\r//gi;}
my @title=split /\t/, $ttt;
while (my $line = <ZZ>){
    chomp $line;
    foreach my $i (keys %searchTerm){
	if (index(lc($line), $i) >=0){
	    my $x=$line;
	    my $j="<font color=red>$i</font>";
	    for ($x) {s/$i/$j/g;}
	    my @ttt=split /\t/, $line;
	    #$searchID{$i}{$ttt[0]}=$j;
	    $search{$ttt[0]}=$i;
	    push(@ID, $ttt[0]);
	}
    }
    my @tmp=split/\t/, $line;
    $ann{$tmp[0]}{"name"}=$tmp[1];
    $ann{$tmp[0]}{"cas"}=$tmp[2];
    $ann{$tmp[0]}{"hmdb"}=$tmp[4];
    $ann{$tmp[0]}{"kegg"}=$tmp[5];
    $ann{$tmp[0]}{"mass"}=$tmp[6];
    $ann{$tmp[0]}{"pubchem"}=$tmp[7];
    $ann{$tmp[0]}{"RI"}=$tmp[10];
    $ann{$tmp[0]}{"path"}=$tmp[12];
    $ann{$tmp[0]}{"subpath"}=$tmp[11];
    
    
   
}
close ZZ;


#######################################
# 2. load stat, and retrieve data
#      link = 'http://iip.pg.com:10080/iip-web/geneSummary.do?&annotationdata=20,56,13,74,44,11,55,50,51,49&acronym=' + gene;

########################################
print "Content-type: text/html\n\n";
print "<HTML><HEAD><TITLE>Metabolite Seach Result</TITLE>
<script language=\"JavaScript\">
</script>\n</HEAD><BODY>\n";

print "<TABLE border=1>\n";
print "<H1>Metabolon Search Result</H1>\n";  
print "<TR><TH>SearchTerm</TH><TH>CompID</TH>";   
foreach my $st (@stats){
    print "<TH>fc $st</TH><TH>pval $st</TH>";
}
print "<TH>ChemName</TH><TH>KEGG</TH><TH>HMDB</TH><TH>PubChemID</TH><TH>CAS</TH><TH>Pathway</TH><TH>SubPathway</TH><TH>MASS</TH><TH>RI</TH></TR>\n";
get_all_stat(\@stats, \@ID, \%search);
print "</TABLE>\n";
print "</BODY></HTML>\n";

sub get_all_stat{
    my ($Hstats, $Hprobe, $Hsearch)=@_;
    my %fc;
    my %pval;
   
    foreach my $st (@{$Hstats}){
	my $stfile=$PATH.$st;
	fill_stat(\%{$fc{$st}}, \%{$pval{$st}}, $stfile);
    }
   
    my %check;
    foreach my $gene (@{$Hprobe}){
	my $kegg=""; if(defined $ann{$gene}{"kegg"}){ $kegg=$ann{$gene}{"kegg"};}
	my $name="";if(defined $ann{$gene}{"name"}){ $name=$ann{$gene}{"name"};}
	my $hmdb="";if(defined $ann{$gene}{"hmdb"}){ $hmdb=$ann{$gene}{"hmdb"};}
	my $pubchem="";if(defined $ann{$gene}{"pubchem"}){ $pubchem=$ann{$gene}{"pubchem"};}
	my $cas="";if(defined $ann{$gene}{"cas"}){ $cas=$ann{$gene}{"cas"};}
	my $path="";if(defined $ann{$gene}{"path"}){ $path=$ann{$gene}{"path"};}
	my $subpath="";if(defined $ann{$gene}{"subpath"}){ $subpath=$ann{$gene}{"subpath"};}
	my $mass="";if(defined $ann{$gene}{"mass"}){ $mass=int($ann{$gene}{"mass"});}
	my $RI="";if(defined $ann{$gene}{"RI"}){ $RI=$ann{$gene}{"RI"};}
	if (!defined $check{$gene}){
	    $check{$gene}=1; #only print once
	    my $searchlink="https://www.google.com/?gws_rd=ssl#q=".$$Hsearch{$gene};
	    my $lline="<TR><TD><a href=\"".$searchlink."\">".$$Hsearch{$gene}."</a></TD><TD><a href=\"$genelink".$gene.".PNG\">".$gene."</a></TD>";
	    my $pp=0;
	    foreach my $st (@{$Hstats}){
		my $FC="NA";
		if(defined $fc{$st}{$gene}){
		    $FC=$fc{$st}{$gene};
		}
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
		my $P="NA";
		if(defined $pval{$st}{$gene}){
		    $P=$pval{$st}{$gene};
		}
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
	    }
	   
	    $lline.="<TH>$name</TH><TH>$kegg</TH><TH>$hmdb</TH><TH>$pubchem</TH><TH>$cas</TH><TH>$path</TH><TH>$subpath</TH><TH>$mass</TH><TH>$RI</TH></TR>";
	    if ($pp ==1 ){print $lline,"\n";}
	}
    }    
    return;
}



####################
# support methods
####################



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
	my @tt=split /\|\|/, $tmp[1];
	for(my $i=1;$i<@tt;$i++){
	    if($tt[$i] ne "---"){
		$$h{$tt[$i]}{$tmp[0]}=1;
		print STDERR $tt[$i],"\t", $tmp[0], "\n"; 
	    }
	}
	$$h{lc($tmp[0])}{$tmp[0]}=1;
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
	if($tmp[2] =~ /\d/){
	    my $ff = sprintf("%.4f", $tmp[2]);
	    $$h_fc{lc($tmp[0])}=$ff;
	}
	if($tmp[1] =~ /\d/){
	    my $pp = sprintf("%.4f", $tmp[1]);
	    $$h_p{lc($tmp[0])}=$pp;
	}
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
