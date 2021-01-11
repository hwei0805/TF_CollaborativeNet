#!/usr/bin/perl
#=============================================================================#
#         FILE:  SCCM_builder.pl
#
#  DESCRIPTION:  read gene id(usually TF gene) from a file, gets its expression vector from the expression file,
#                Computes the rank correlation with each vector(line) in expression file. 
#                For each TF gene choose top 100 most correlated genes,builds a TF*TF matrix the value is the shared top 100 genes.   
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  The title line of both files should start with #
#       AUTHOR:  Jeff Nie, <jnie@morgridgeinstitute.org>>
#      COMPANY:  Morgridge Research Institute
#      VERSION:  0.0.1
#      CREATED:  10/06/10 11:14:40 CST
#     REVISION:  ---
#===============================================================================
use Statistics::RankCorrelation;
use File::Basename;
use Parallel::ForkManager;
use strict;
my %args=@ARGV;
my $genelist = $args{'-g'};
my $expression = $args{'-e'};
my $top=$args{'-t'}||100;
my $cpu=$args{'-cpu'}||1;
unless($genelist and $expression){
    print "Usage: 
    $0 -g genelist_file -e expression_file <-t NumOfTopGeneChoose> <-cpu howManyCPUs> \n";
    exit;
}
my $expbasename = basename($expression);
my $out=$genelist."_".$expbasename."_coexp.matrix.txt";
open (L,"$genelist")|| die "can't open genelist file $genelist $!";
open (EXP,"$expression")|| die "can't open expression file $expression $!";
open (OUT,">$out")|| die "can't open output file $out $!";
print "SCCM pipeline start from ",`date`;
my %hash;
my @id_list;
my %tf;
my @TFgene;
while(<EXP>){
    chomp;
    next if /^\#/;
    my @array = split;
    my $id = shift @array;
    push (@id_list,$id);
    $hash{$id}= [@array];
}
my $pm=new Parallel::ForkManager($cpu);
my $test=0;
mkdir "top_$top", 0777 unless -d "top_$top";

while(<L>){
    chomp;
    my $gene =$_;
    push @TFgene,$gene;
    $pm->start and next;
    next if /^\#/;
    my %rho;
    foreach (@id_list){
        my $c = Statistics::RankCorrelation->new( $hash{$gene}, $hash{$_} );
	my $n = $c->spearman;
	$n=sprintf("%.4f",$n);
	$rho{$_}=$n;
    }
    #print `date`,"\n"; 
    my @sorted = sort {$rho{$b} <=> $rho{$a}} keys %rho; #list gene names in rho sorted order
    open (TOP,">top_$top/$gene")|| die "can't open output file top_$top/$gene $!";
    for(1..$top){
       print TOP "$sorted[$_]\n";
    }     
    close TOP;
    $pm->finish;
}
$pm->wait_all_children;
close L;
print "Correlation done at ",`date`,"\n" ;

opendir(DIR, "top_$top") or die "can't opendir top_$top: $!";
while (defined(my $file = readdir(DIR))) {
    next if $file =~ /^\.\.?$/; # skip . and ..
    #push @TFgene, $file;
    open (TF,"top_$top/$file")|| die "can't open input file $file $!";	
    while(<TF>){
	chomp;
        $tf{$file}{$_}=1; 
    }
    close TF;
}
closedir(DIR);
print "Read top correlated TF gene files done at ",`date`;
foreach my $rGene (@TFgene){ # will be the row TF gene
    print  OUT "$rGene\t";
    foreach my $cGene (@TFgene){ # will be the col TF gene
	my $count=0;
	foreach my $k (keys %{$tf{$cGene}}){ #compare top 100 genes of col TF gene with top 100 row TF genes see how many of them are the same
	   $count ++ if $tf{$rGene}{$k};
	}
	print OUT "$count\t";
    }
    print OUT "\n";
}

print "The program end at: ",`date`;


=comment start
# ########## format ##############
genelist file: 
the id  should have same format with geneid in expression file
and all id in genelist should be included in expression file

#geneid
NM_012138&AATF
BC098116&ABCA11
BC105084&ABCA11
BC099730&ABCA11
NM_013375&ABT1
NM_001616&ACVR2A
NM_015339&ADNP
NM_031284&ADPGK

#expressioin file: tab delimited file

#geneid          day1   day2    day3   day4  day5   day6   day8  day9   day10  day11  day23   day13
B000409&MKNK1    0.685  0.304  0.391  0.564  0.352  0.463  0.64  0.291  0.535  0.373  0.588  0.6344
AB000461&C4orf8  0.148  0.591  1.011  0.398  0.784  1.315  0.69  1.753  0.671  0.879  0.702  1.392

# output format

AB000409&MKNK1      BX640837&SLC26A3  0.3262
AB000461&C4orf8     BX640837&SLC26A30  0.0733
AB000781&KNDC1      BX640837&SLC26A300  0.0676
AB001025&RYR3       BX640837&SLC26A3000  0.1418

=cut

