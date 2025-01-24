#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH --array=0-15
#SBATCH --output=resAlignSTAR.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=6000

# use 4 cores for each node, as many nodes as samples to align
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

# load STAR module
ml releases/2019b STAR/2.7.3a-GCC-8.3.0

# set genome index directory
indexDir=/scratch/ulb/BCTL/apapagia/resources/refgenie/data/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index/default/

# find fastq files
fastqDir=/scratch/ulb/BCTL/apapagia/neoCell
fq1=($(find $fastqDir -type f -name "*_R1_001.fastq.gz" | sort))
fq2=($(find $fastqDir -type f -name "*_R2_001.fastq.gz" | sort))
# retrieve sample names
id=$(basename ${fq1[$SLURM_ARRAY_TASK_ID]} _R1_001.fastq.gz)

# store data in local scratch creating a directory for each sample/node
SCRATCH=$LOCALSCRATCH/$SLURM_JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

# run STAR align
printf "Processing reads\nR1: ${fq1[$SLURM_ARRAY_TASK_ID]} \nR2: ${fq2[$SLURM_ARRAY_TASK_ID]} \n\n ID: $id"
STAR --outFileNamePrefix $SCRATCH/$id --outSAMtype BAM SortedByCoordinate --genomeDir $indexDir --runThreadN 8 \
	--readFilesCommand zcat --readFilesIn ${fq1[$SLURM_ARRAY_TASK_ID]} ${fq2[$SLURM_ARRAY_TASK_ID]} \
	--twopassMode Basic

# copy data to desired output directory
outDir=/scratch/ulb/BCTL/apapagia/star_align
mkdir $outDir
srun cp -r  $SCRATCH/*  $outDir || exit $?
srun rm -rf $SCRATCH || exit $?