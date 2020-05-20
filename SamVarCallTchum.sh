#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=sam_varcall
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu


#SAMTOOLs mpileup version 1.5 options used:
#C = adjust mapping quality; recommended:50, disable:0 [0]
#d = max per-file depth; avoids excessive memory usage [250]
#f = faidx indexed reference sequence file
#q = skip alignments with mapQ smaller than INT [0]
#Q = skip bases with baseQ/BAQ smaller than INT [13]
#g = generate genotype likelihoods in BCF format
#t = --output-tags LIST  optional tags to output:DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR []

#BCFTOOLs call version 1.6 options used
#v = output variant sites only
#c/m = he original calling method (conflicts with -m) or alternative model for multiallelic and rare-variant calling (conflicts with -c)
#p = variant if P(ref|D)<FLOAT with -c [0.5]
#P =  --prior <float-o mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]
#O =  output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' (here it is 'v') 
#o = write output to a file [standard output]

module load samtools
module load bcftools

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/timema_chumash_gbs/Alignments_mel

samtools mpileup -C 50 -d 250 -f /uufs/chpc.utah.edu/common/home/u6000989/data/timema/tcrDovetail/version3/mod_map_timema_06Jun2016_RvNkF702.fasta -q 20 -Q 20 -g -I -t DP,AD -u -b bamfiles -o tchum_mel_gbs_2019.bcf

bcftools call -v -c -p 0.01 -P 0.001 -O v -o tchum_mel_gbs_2019.vcf tchum_mel_gbs_2019.bcf




