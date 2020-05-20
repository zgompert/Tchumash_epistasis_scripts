#!/usr/bin/perl
#
# run gemma jobs
#

use Parallel::ForkManager;
my $max = 20;
my $pm = Parallel::ForkManager->new($max);
my $ph = "pheno_color.txt";

foreach $g (@ARGV){

	$base = $g;
	$base =~ s/mod_g_//;
	$base =~ s/\.txt//;
	foreach $ch (0..9){
		system "sleep 2\n";
		foreach $trait (1..2){
			$pm->start and next; ## fork
			$out = "oe_$base"."_ph$trait"."_ch$ch";
			system "gemma -g $g -p $ph -bslmm 1 -n $trait -k output/tchum_epi.cXX.txt -notsnp -w 200000 -s 1000000 -o $out\n";
			$pm->finish;
		}
	}
}

$pm->wait_all_children;



