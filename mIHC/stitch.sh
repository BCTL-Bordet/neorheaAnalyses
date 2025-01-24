#!/bin/bash
#SBATCH --job-name=qupath
#SBATCH --array=0-14
#SBATCH --output=logs/qupath/stitch_%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=26000

# set path to QuPath executable, script to run and directories of samples
quPathPath=~/software/QuPath/bin/QuPath
scriptToRun=~/qupathStitch.groovy
dataDirs=~/data/for_upload_3_reload/Monoplexes_for_upload_in_ceci/

# find directories where unstitched tifs are located
dirs=($(find $dataDirs -maxdepth 1 -mindepth 1	 -type d))

# find first tif file for each of the directories
smpls=(); for dir in ${dirs[@]}; do smpls+=($dir"/"$(ls $dir | head -n 1)); done

# run qupath stitch script for each of the samples
$quPathPath script $scriptToRun -i ${smpls[$SLURM_ARRAY_TASK_ID]}