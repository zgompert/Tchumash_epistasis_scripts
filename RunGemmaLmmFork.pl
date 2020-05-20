#!/usr/bin/perl
#
# run gemma jobs
#

use Parallel::ForkManager;
my $max = 12;
my $pm = Parallel::ForkManager->new($max);


foreach $g (@ARGV){
	$ph = $g;
	$ph =~ s/mod_g/ph/;
	$base = $g;
	$base =~ s/mod_geno_//;
	$base =~ s/\.txt//;
	foreach $trait (1..3){
		$pm->start and next; ## fork
		$out1 = "k_$base"."_ph$trait";
		$out2 = "o_$base"."_ph$trait";
		system "gemma -g $g -p $ph -gk -o $out1\n";
		system "gemma -g $g -p $ph -lmm 4 -k output/$out1".".cXX.txt -n $trait -maf 0.0 -o $out2\n";
		$pm->finish;
	}
}

$pm->wait_all_children;



