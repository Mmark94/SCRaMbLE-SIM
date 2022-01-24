#! /usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Cwd;

my $usage=<<USAGE;
Describe:Only used to refactor SCRaMbLEd genome sequence.
  
  Author: Wang Yun, wangyun\@genomics.org.cn

 Version: 1.0.1  Date: 2021-8-19
    Note: only get all sequence of the path in <encod> file.

   Usage:perl $0 <ref_fa> <loxpregion> <encod> [option]
         
         -spid     <str>  sample id
         -outdir          output directory
         -help            show this message
         
 History Version:
     1.0  Date: 2012-11-05 # only get one sequence of the first path in <encod> file.

USAGE
  
my ($sample,$reseqid,$outdir,$help);
GetOptions(
   "spid:s" => \$sample,
   "outdir:s" => \$outdir,
   "help" => \$help,
);

die $usage if (@ARGV != 3 ||$help);

my ($ref_fa,$loxpregion,$encod) = @ARGV;

open IN1, $ref_fa ||die $!;
open IN2, $loxpregion ||die $!;
open IN3, $encod || die $!;

$outdir ||= getcwd;
open OUT,">$outdir/$sample.reseq.fasta" || die $!;

#print "loading data...\n";

my %lxpregseq;
my $refseq=();
my $n=0;
my %_idcheck;

<IN1>;
while (<IN1>){
	chomp;
    $_ =~ s/\r//ig;
	$refseq .= $_;
}

my $len= length($refseq);
print "Reference length: $len\n";

while(<IN2>){
	chomp;
	my ($star,$END,$indx)=(split /\t/)[3,4,8];
	my $reglen = $END - $star +1;
	 
	$lxpregseq{$indx} = substr($refseq,$star-1,$reglen);
    $n++;
 }
print "Total $n loxp region!\n";	

while (<IN3>){
   chomp $_;
   my $restructseq = $_;

   chomp $restructseq;
   my ($solution_id, $encode) = (split /\s+/,$restructseq)[0,1];
   my @codeord = split /,/,$encode;

   my $reseq;
   foreach (@codeord){
	   if ($_ == 0){
		   next;
	    }elsif($_>0){
	       $reseq .= $lxpregseq{$_};
	    }elsif($_<0){
	       my $idx = abs($_);
	       my $revcom = reverse($lxpregseq{$idx});
	       $revcom =~ tr/ATGCatgc/TACGtacg/;
	 	   $reseq .= $revcom;
  	    }
   }

   print OUT ">$solution_id\n";
   print OUT "$reseq\n";

   my $lenscb = length($reseq);
   print "$solution_id Refactoring sequence length:$lenscb\n";
}

close IN1;
close IN2;
close IN3;
close OUT;

#===end===	
