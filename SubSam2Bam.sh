#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=samtools
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/timema_chumash_gbs/Alignments_mel

perl Sam2BamFork.pl *sam

