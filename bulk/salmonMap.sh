#!/bin/bash
#SBATCH --job-name=salmonMap
#SBATCH --array=0-15
#SBATCH --output=resSalmonMap.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=8000

# set path of salmon
salmonDir=/home/ulb/BCTL/apapagia/software/salmon-1.5.2_linux_x86_64/bin
# set path of salmon genome index
salmonIdx=/scratch/ulb/BCTL/apapagia/resources/refgenie/data/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/salmon_sa_index/default/

# find paired forward fastq read files
fastqDir=/scratch/ulb/BCTL/apapagia/trimmomatic_out_neoRhea
fq=($(find $fastqDir -type f \( -name "*1P.fq.gz" -o -name "*_1P.fastq.gz.fq.gz" \)))

# get ID of sample
id=$(basename ${fq[$SLURM_ARRAY_TASK_ID]} _1P.fq.gz)
printf "Processing read\nR1: ${fq[$SLURM_ARRAY_TASK_ID]} \n"

# store data in local scratch creating a directory for each sample/node
SCRATCH=$LOCALSCRATCH/$SLURM_JOB_ID
mkdir -p $SCRATCH

# run salmon single-end mode
cd $salmonDir
./salmon quant -i $salmonIdx -l A \
    -r ${fq[$SLURM_ARRAY_TASK_ID]} \
    -o $SCRATCH/$id \
    --validateMappings --gcBias --threads 4

# copy data to desired output directory
outDir=/scratch/ulb/BCTL/apapagia/salmon_MCF7
mkdir $outDir
srun cp -r  $SCRATCH/*  $outDir || exit $?
srun rm -rf $SCRATCH || exit $?