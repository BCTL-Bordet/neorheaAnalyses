#!/bin/bash
#SBATCH --job-name=salmonIdx
#SBATCH --output=resSalmonIdx.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=8000

# set salmon directory
salmonDir=/home/ulb/BCTL/apapagia/software/salmon-1.5.2_linux_x86_64/bin
# set path of transcripts.fa file
hg38fa=/scratch/ulb/BCTL/apapagia/resources/gencode/gencode.v38.lncRNA_transcripts.fa
# set desired output directory
outDir=/scratch/ulb/BCTL/apapagia/resources/salmon_hg38Idx_lncRNAs

cd $salmonDir
./salmon index -t $hg38fa -i $outDir --gencode