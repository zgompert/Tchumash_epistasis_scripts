#!/usr/bin/perl
#
# run gemma jobs
#

use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);
my $ph;
my @phf = ("pheno_rg_cv.txt","pheno_gb_cv.txt");
my $trt;

foreach $g (@ARGV){
	foreach $ph (@phf){
		$base = $g;
		$base =~ s/mod_g_//;
		$base =~ s/\.txt//;
		$ph =~ m/pheno_([a-z]+)/;
		$trt = $1;
		foreach $ch (0..4){
			system "sleep 2\n";
			foreach $set (1..10){
				$pm->start and next; ## fork
				$out = "o_cv_$base"."_ph$trt"."_set$set"."_ch$ch";
				if($base eq "tchum"){
					system "gemma -g $g -p $ph -bslmm 1 -n $set -maf 0.0 -w 200000 -s 1000000 -o $out\n";
				}
				elsif($base eq "epi_tchum"){
					system "gemma -g $g -p $ph -bslmm 1 -n $set -k output/tchum_epi.cXX.txt -notsnp -maf 0.0 -w 200000 -s 1000000 -o $out\n";
				}
				$pm->finish;
			}
		}
	}
}

$pm->wait_all_children;



