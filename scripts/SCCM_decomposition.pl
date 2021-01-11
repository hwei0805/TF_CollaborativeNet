#!/usr/bin/perl
#===============================================================================
#
#         FILE:  SCCM_decomposition.pl
#
#        USAGE:  SCCM_decomposition.pl -g genelist_file -m matrix_file <-ks KickOutSize> <-tl1 linkCutOff1> <-tl2 linkCutOff2> <-tl3 linkCutOff3>
#
#  DESCRIPTION:  This program reads a file of gene name list and a
#                matrix of correlation between genes, arranges the
#                genes into groups and outputs a file of list of
#                the groups
#
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Yao Xue        08/05/2010
#      VERSION:  0.0.1
#      CREATED:  09/05/2010 11:14:40 CST
#     Modified:  1. Jeff Nie, 10/01/2010 add serveral parameters and clean up some code. 
#     REVISION:  ---
#===============================================================================


use File::Basename;
use strict;
my %args = @ARGV;
my $genelist = $args{'-g'};
my $matrix = $args{'-m'};
my $kickedSize=$args{'-ks'}||2000;
my $tl1=$args{'-tl1'}||"1.5";
my $tl2=$args{'-tl2'}||"1.2";
my $tl3=$args{'-tl3'}||"0.8";
my $i;
my $j;
my @kickedOut;
my @genes;
my @finList;
my @sorted;
unless($genelist and $matrix){
    print "Usage:
    $0 -g genelist_file -m matrix_file <-ks KickOutSize> <-tl1 linkCutOff1> <-tl2 linkCutOff2> <-tl3 linkCutOff3> \n";
    exit;
}
print "$0 start at ", `date`;
my $expbasename = basename($matrix);
my $out="Cluster_".$expbasename;
open (L,"$genelist")|| die "can't open genelist file $genelist $!";
open (EXP,"$matrix")|| die "can't open matrix file $matrix $!";
open (OUT,">$out")|| die "can't open output file $out $!";
my (%hash, %hashcp, @inLink, $dir, $tmp, $inlen);
my @id_list;
my $listsize;
my $sqrsum = 0;
my @calc;
my $avg;
my $std;
my $sum = 0;
my $clen;
$dir = "cluster_within_connectivity";
if (-d "./$dir") {
    qx(rm -r $dir);
}
qx(mkdir $dir);

$i=0;
while(<L>){
    chomp();
    $genes[$i] = $_;
    print "gene: $i\t$genes[$i]\n";
    $i++;
}

$listsize = @genes; 
while(<EXP>){
    chomp;
    next if /^\#/;
    my @array = split(/\t/, $_);
    my $id = shift @array;
    for($i=0;$i<$listsize;$i++){
	if($id eq $genes[$i]){
	    $hash{$id.$genes[$i]}=0;
	} else {
	    $hash{$id.$genes[$i]}= $array[$i];
            $hashcp{$id.$genes[$i]}= $array[$i];  # keep a  copy of non-zero TF pair for later use 
	    if($array[$i]>0){
		push(@calc, $array[$i]);
	    }

	}
    }
}
print "Matrix file processing done at ",`date`;
foreach(@calc){
    $sum += $_;
}
$clen=@calc;
print "clen: $clen\n";
$avg = $sum/$clen;

foreach(@calc){

    $sqrsum += ($_ - $avg)**2;
}
$std = ($sqrsum/@calc)**0.5;

printf ("avg = $avg\n");
printf ("std = $std\n");

my $lv1 = $avg + $tl1*$std;
my $lv2 = $avg + $tl2*$std;
my $lv3 = $avg + $tl3*$std;
my $lvb = $std/2;

@sorted = sort {$hash{$b} <=> $hash{$a}} keys %hash;  # sort to get the pair of TFs with maximal value
print "sort0:$sorted[0]\n";

#make a subdir to store the results

my $nc=1;
push (@finList, "cluster: $nc \n");
while($hash{$sorted[0]} > $lv3){
    print "sort_2nd_0:$sorted[0]\n";
 
    $listsize = @genes;
     
     # search the pair of TFs in @genes and delete them while adding them to the @kickedOut
    ### why not directly deleted from %hash, for delete from @gene?(we could do foreach @gene set i and j gene to '')
    for($i=0;$i<$listsize;$i++){ # this whole loop is just for the purpose get ride of the kickout pair from hash 
	for($j=0;$j<$listsize;$j++){
	    if(($genes[$i].$genes[$j] eq $sorted[0]) || ($genes[$j].$genes[$i] eq $sorted[0])){
		push(@kickedOut,$genes[$i]);
		push(@kickedOut,$genes[$j]);
		delete $hash{($genes[$i].$genes[$j])};
		delete $hash{($genes[$j].$genes[$i])};
		$genes[$i] = '';   # set $genes[i] to 0 though not deleted
		$genes[$j] = '';   # set $gene [j] = '' though not deleted
	    }
	}
    }
    print "kickout: @kickedOut\n";    # print the pair of TFs
    while ($kickedSize != @kickedOut){
	$kickedSize = @kickedOut;
	printf "$kickedSize\n";
	$listsize = @genes;
	for($i=0;$i<$listsize;$i++){
	    my $high;
	    my $mid;
	    my $low;
	    if($kickedSize < 3){  # make they third TF enter to cluster easier
		$high = 0;
		$mid = 0;
		$low = 1;
	    }
	    else{
		$high = 0;
		$mid = 0;
		$low = 0;
	    }
	    foreach(@kickedOut){
		if($hash{($_.$genes[$i])} < $lv3 ){
		    last;
		}
		if(($hash{($_.$genes[$i])} > $lv1) && ($high == 0)){
		    $high = 1;
		}
		elsif(($hash{($_.$genes[$i])} > $lv2) && ($mid == 0)){
		    $mid = 1;
		}
		elsif(($hash{($_.$genes[$i])} > $lv3) && ($low == 0)){
		    $low = 1;
		}
		if($high && $mid && $low ){
		    push (@kickedOut,$genes[$i]);   # add a linked gene
                    printf "$genes[$i]\n";
		    $genes[$i] = '';   # set this gene name as '' but not delete the slot
		    last;
		}

	    }
	}
    }
    #by now the one cluster building is done
   # now copy these genes in @kickedOut to @finList, while deleting all of them from %hash 
    #push (@finList, "cluster: $nc \n");
    my $initSize = @kickedOut;

    # before deletion, figure out within connectivity of a cluster
    for($i=0;$i<$initSize;$i++){
	for($j=0;$j<$initSize;$j++){
            if (($hashcp{$kickedOut[$i].$kickedOut[$j]} > $lv3) || ($hashcp{$kickedOut[$j].$kickedOut[$i]} > $lv3)) {
                 push @inLink, {TF1 => $kickedOut[$i], TF2 => $kickedOut[$j], Share => $hashcp{$kickedOut[$i].$kickedOut[$j]} };
		 push @inLink, {TF1 => $kickedOut[$j], TF2 => $kickedOut[$i], Share => $hashcp{$kickedOut[$j].$kickedOut[$i]} };
            }    
	}
    }
    # save clusterlink
    $tmp = "cluster_" . $nc . "_within_link.txt";
    print "tmp: $tmp\n";
    open (CL,">./$dir/$tmp")|| die "can't open output file $out $!";
    $inlen=@inLink;
    for ($i=0; $i<$inlen;$i++){ #print the cluster connectivity file
	print CL "$inLink[$i]{'TF1'}\t$inLink[$i]{'TF2'}\t$inLink[$i]{'Share'}\n";
    }   
    close(CL);

    # delete all paird genes: one is a TF from kickedOut the other is one from @genes and @finList
    for($i=0;$i<$initSize;$i++){
	foreach(@genes){     # for any genes that paired with TFs in kickedOut, delete them.
	    delete $hash{($kickedOut[0].$_)};  
	    delete $hash{($_.$kickedOut[0])};
	}
	foreach(@finList){  # for any TFs that paired with TFs in finList, delete them:  this may be redundent
	    delete $hash{($kickedOut[0].$_)}; 
	    delete $hash{($_.$kickedOut[$0])};
	}
	my $temp = shift @kickedOut;  # store in finList
	push (@finList, $temp);
	$kickedSize = @kickedOut;   # update the size
    }
 
    undef(@kickedOut);   # free KickedOut

    $nc++;
    push (@finList,"cluster: $nc");
    printf "cluster: $nc \n";
    @sorted = sort {$hash{$b} <=> $hash{$a}} keys %hash;   # sort hash again
}


foreach(@finList){ # print the cluster file
    print OUT "$_\n";
  
}

print "$0 done at ",`date`;
close(L);
close(EXP);
close(OUT);
