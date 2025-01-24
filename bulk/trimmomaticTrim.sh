#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --array=0-15
#SBATCH --output=resTrimmomatic.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=8000

# use 4 cores for each node
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

# find fastq files
fastqDir=/scratch/ulb/BCTL/apapagia/MCF7/210309_A00154_0728_BHYC3VDMXX_MIG_u2-hs_RNAseq/
fq1=($(find $fastqDir -type f \( -name "*_R1_001.fastq.gz" \) | sort))
fq2=($(find $fastqDir -type f \( -name "*_R2.fastq.gz" -o -name "*_R2_001.fastq.gz" \) | sort))

# retrieve sample name
id=$(basename ${fq1[$SLURM_ARRAY_TASK_ID]} _R1_001.fastq.gz)

# store data in local scratch creating a directory for each sample/node
SCRATCH=$LOCALSCRATCH/$SLURM_JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

# run trimmomatic
printf "Processing reads:\nR1: ${fq1[$SLURM_ARRAY_TASK_ID]} \nR2: ${fq2[$SLURM_ARRAY_TASK_ID]} \n\n"
java -jar /home/ulb/BCTL/apapagia/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
    -threads 4 -trimlog $SCRATCH/$id"_log.txt" \
    ${fq1[$SLURM_ARRAY_TASK_ID]} ${fq2[$SLURM_ARRAY_TASK_ID]} \
    -baseout $SCRATCH/$id.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:1

# copy data to desired output directory
outDir=/scratch/ulb/BCTL/apapagia/trimmomatic_out_neoRhea
mkdir $outDir
srun cp -r  $SCRATCH/*  $outDir || exit $?
srun rm -rf $SCRATCH || exit $?