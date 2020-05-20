#!/usr/bin/perl

use warnings;
use strict;

# filter vcf files

# usage vcfFilter.pl infile.vcf
#
#
#### stringency variables, edits as desired
# T. chumash N = 454
my $minCoverage = 908; # minimum number of sequences; DP, no. of individuals = X; mincoverage = 2X
my $minAltRds = 10; # minimum number of sequences with the alternative allele; AC
my $notFixed = 1.0; # removes loci fixed for alt; AF
my $mq = 30; # minimum mapping quality; MQ
my $miss = 181; # maximum number of individuals with no data = 20%
my $minMaf = 4; # minimum MAF, as allele count ~ 0.005
my $ratDifDir = 0.01; ## maximum proportion of reads in reverse (relative to avg.) orientation
my $d;
my $frd;
my $tot;

my @line;

my $in = shift(@ARGV);
open (IN, $in) or die "Could not read the infile = $in\n";
$in =~ m/^([a-zA-Z0-9_]+\.vcf)$/ or die "Failed to match the variant file\n";
open (OUT, "> filtered2x_$1") or die "Could not write the outfile\n";

my $flag = 0;
my $cnt = 0;

while (<IN>){
	chomp;
	$flag = 1;
	#print "\n";
	if (m/^\#/){ ## header row, always write
		$flag = 1;
	}
	elsif (m/^ScRvNkF_/){ ## this is a sequence line, you migh need to edit this reg. expr.
		$flag = 1;
		$d = () = (m/:0,0,0:0:0,0/g);
		if ($d >= $miss){
			$flag = 0;
	#		print "fail missing : ";
		}
		else{
	#		print "pass missing : ";
		}
		if (m/[ACTGN]\,[ACTGN]/){ ## two alternative alleles identified
			$flag = 0;
	#		print "fail allele : ";
		}
		@line = split(/\s+/,$_);
		if(length($line[3]) > 1 or length($line[4]) > 1){
			$flag = 0;
	#		print "fail INDEL : ";
		}
		m/DP=(\d+)/ or die "Syntax error, DP not found\n";
		if ($1 < $minCoverage){
			$flag = 0;
	#		print "fail DP : ";
		}
		m/DP4=\d+,\d+,(\d+),(\d+)/ or die "Syntax error, DP4 not found\n";
		if (($1+$2) < $minAltRds){
			$flag = 0;
	#		print "fail alt : ";
		}
		m/DP4=(\d+),(\d+),(\d+),(\d+)/ or die "Syntax error, DP4 not found\n";
		$frd = $1 + $3;
		$tot = $1 + $2 + $3 + $4;
		if (($frd/$tot) > $ratDifDir and ($frd/$tot) < (1-$ratDifDir)){
			$flag = 0;
	#		print "fail orientation : ";
		}
		if(m/AF1=1;/){
			$flag = 0;
		}
		#my file has AC1 instead of AC, so I have changed this - Sam	
		m/AC1=(\d+)/ or die "Syntax error, AC not found\n";
		if ($1 < $minMaf or $1 > (2840 - $minMaf)){
			$flag = 0;
	#		print "fail AC : ";
		}
		m/MQ=([0-9\.]+)/ or die "Syntax error, MQ not found\n";
		if ($1 < $mq){
			$flag = 0;
	#		print "fail MQ : ";
		}
		if ($flag == 1){
			$cnt++; ## this is a good SNV
		}
	}
	else{
		print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
		$flag = 0;
	}
	if ($flag == 1){
		print OUT "$_\n";
	}
}
close (IN);
close (OUT);

print "Finished filtering $in\nRetained $cnt variable loci\n";
