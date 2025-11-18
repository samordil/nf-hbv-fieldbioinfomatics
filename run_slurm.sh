#!/bin/bash
#SBATCH --job-name=hbv
#SBATCH --partition=ncpu
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=07-00:00:00
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@email.com

export NXF_WORK=/scratch/$USER/nxf_work
export TMPDIR=/scratch/$USER/tmp
mkdir -p "$NXF_WORK" "$TMPDIR"

# Activate existing conda environment
source ~/Tools/miniforge3/etc/profile.d/conda.sh
conda activate hbv-artic

# Run the pipeline
nextflow run main.nf -profile hbv,slurm -resume
