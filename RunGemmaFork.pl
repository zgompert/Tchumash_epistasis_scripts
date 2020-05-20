#!/usr/bin/perl
#
# run gemma jobs
#

use Parallel::ForkManager;
my $max = 20;
my $pm = Parallel::ForkManager->new($max);


foreach $g (@ARGV){
	$ph = $g;
	$ph =~ s/mod_g/ph/;
	#$ph =~ s/\.txt/_dw.txt/;
	$base = $g;
	$base =~ s/mod_geno_//;
	$base =~ s/\.txt//;
	foreach $ch (0..9){
		system "sleep 2\n";
		#foreach $trait (1..2){
		foreach $trait (1..3){
			$pm->start and next; ## fork
			#$out = "o_$base"."_dw_ph$trait"."_ch$ch";
			#$out = "o_$base"."_ph$trait"."_ch$ch";
			$out = "o_$base"."_ph$trait"."_ch$ch";
			system "gemma -g $g -p $ph -bslmm 1 -n $trait -maf 0.0 -w 200000 -s 1000000 -o $out\n";
			$pm->finish;
		}
	}
}

$pm->wait_all_children;



