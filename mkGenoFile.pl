#!/usr/bin/perl
#
# formats the geno file
#

## get actual SNP ids
open(IN, "tpetSnps.txt") or die "failed to read SNPs file\n";
#open(IN, "Snps.txt") or die "failed to read SNPs file\n";
while(<IN>){
	chomp;
	push(@snps,$_);
}
close(IN);

foreach $in (@ARGV){ ## genoe files

	## read and write geno file
	open(IN, $in) or die "failed to read the genotype file\n";
	$out = "mod_$in";
	open(OUT, "> $out") or die "failed to write\n";
	$i = 0;
	while(<IN>){
		chomp;
		s/ /, /g or die "failed space sub\n";
		print OUT "$snps[$i], A, T, $_\n";
		$i++;
	}
	close(IN);
	close(OUT);
}
