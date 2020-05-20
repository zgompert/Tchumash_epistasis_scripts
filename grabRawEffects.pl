#!/usr/bin/perl

## read in snps and create hash
open(IN,"snpList.txt") or die "failed to read the snp file\n";
while(<IN>){
	chomp;
	push(@snps,$_);
	$eff{$_} = 0;
}
close(IN);

print "finished reading snp list\n";

## combine information over reps and sort
$in = shift(@ARGV); ## 0 rep file, from fork script
open(OUT, "> rawbeta_$in") or die "failed to write mav $in\n";
$nreps=10;
foreach $rep (0..9){
	$in =~ s/ch\d+/ch$rep/ or die "failed sub for $in to $rep\n";
	print "working on $in\n";
	open(IN, $in) or die "failed to open the infile $in\n";
	<IN>; ## burn header
	while(<IN>){
		chomp;
		@line = split(/\s+/,$_);
		$eff{$line[1]} += ($line[5]) * 0.1;## * .1 to divide by 10
	}
	close(IN);
}

foreach $snp (@snps){
	print OUT "$eff{$snp}\n";
}
close(OUT);
