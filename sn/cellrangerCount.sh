#!/bin/bash
#SBATCH --job-name=cRangerC
#SBATCH --array=0-66
#SBATCH --output=logs/cellranger-arc_count_%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=6000

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

hg38ref=/scratch/ulb/BCTL/apapagia/resources/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

libCSVsDir=/scratch/ulb/BCTL/apapagia/sc_RNA_ATAC/cellranger/cellrangerLibraryCSVs
libCSVs=($(find $libCSVsDir -type f -name "*.csv"))

sampleFile=${libCSVs[$SLURM_ARRAY_TASK_ID]}
id=$(basename $sampleFile .csv)

# be a good boi and store data in local scratch of each node
SCRATCH=$LOCALSCRATCH/$SLURM_JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH
outDir=$SCRATCH


# run only if previous iteration didn't result into anythin
if [ ! -d $outDir/$id/outs ]; then
        cellranger-arc count --id=$id \
                --reference=$hg38ref \
                --libraries=$sampleFile \
                --no-bam \
                --disable-ui \
                --localcores=8 \
                --localmem=46
fi

# copy data to global scratch
srun cp -r  $SCRATCH/* /scratch/ulb/BCTL/apapagia/sc_RNA_ATAC/cellranger/cellrangerOut/  || exit $?
srun rm -rf $SCRATCH || exit $?