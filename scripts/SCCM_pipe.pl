#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  SCCM_Pipe.pl
#
#        USAGE:  ./SCCM_Pipe.pl <config_file>
#
#  DESCRIPTION:  pipeline to run the SCCM package
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Jeff Nie, <jnie@morgridgeinstitute.org>>
#      COMPANY:  Morgridge Research Institute
#      VERSION:  0.0.1
#      CREATED:  10/02/10 11:14:40 CST
#     REVISION:  ---
#===============================================================================
use strict;
use File::Basename;
use Config::Simple;
if($#ARGV < 0){
   print "Usage:
   $0 path_of_config_file
   example:
   $0 /home/jnie/ES.cfg\n";
   exit;
}

print  "$0 @ARGV\n";
print "Pipeline  started at ",`date`,"...\n";
qx(echo $0 @ARGV >standout);
qx(echo Pipeline  started at >>standout);
qx(date >>standout);

my $cfg = new Config::Simple(shift(@ARGV)) or die Config::Simple->error();
my $lib=$cfg->param('lib');# directory that hold all programs called by this programs. 
my $cpu=$cfg->param('cpu');# how many CPUs you want to use. 
my $genelist=$cfg->param('geneList');
my $expression=$cfg->param('expression');
my $expbasename = basename($expression);
my $top=$cfg->param('topPick');	
#my $skipBuilder=$cfg->param('skipBuilder');
my $kickedSize=$cfg->param('kickSize');
my $tl1=$cfg->param('tripleLink1');
my $tl2=$cfg->param('tripleLink2');
my $tl3=$cfg->param('tripleLink3');
my $indir=dirname($genelist);
my $outdir=$cfg->param('outdir')||$indir;
my $matrix=$cfg->param('matrixName')||$genelist."_".$expbasename."_coexp.matrix.txt";
open(LOG, ">$outdir\/log") ||die "can't open log file. $!";
print  LOG "$0 @ARGV\n";
print LOG "Pipeline  started at ",`date`,"...\n";
qx(echo start building co-expression matrix ... >>standout);
qx($lib/SCCM_builder.pl -g $genelist -e $expression -t $top -cpu $cpu >>standout) unless (-s $matrix);#skip SCCM building step if the matrix file already exist.
qx(echo end building co-expression matrix. >>standout);
qx(echo start matrix decomposition ... >>standout);
print LOG "SCCM_builder.pl finished  at ",`date`,"...\n";
qx($lib/SCCM_decomposition.pl -g $genelist -m $matrix -ks $kickedSize  -tl1 $tl1 -tl2 $tl2 -tl3 $tl3 >>standout);
print "$lib/SCCM_decomposition.pl -g $genelist -m $matrix -ks $kickedSize  -tl1 $tl1 -tl2 $tl2 -tl3 $tl3\n";
print "Program  finised at ",`date`,"\n";
print LOG "Pipeline finised at ",`date`,"\n";
qx(echo end matrix decomposition >>standout);
qx(echo Pipeline  finished at >>standout);
qx(date >>standout);
my ($user,$systems,$cuser,$csystem)=times();
printf LOG "This pipeline has consumed tolal %.3f seconds\n", $user+$systems+$cuser+$csystem;
close LOG;
