#!/bin/bash
#SBATCH --job-name=hbv
#SBATCH --partition=highmem
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=07-00:00:00
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=use_your_email@gmail.com    # Change to your email

export NXF_WORK=/scratch/hbv_nxf_work      # Change to your instituion scratch dir
export TMPDIR=/scratch/hbv_tmp             # Change to your instituion scratch dir
mkdir -p "$NXF_WORK" "$TMPDIR"

# Activate existing conda environment
source ~/miniforge3/etc/profile.d/conda.sh   # Change to point ot you custom environment
conda activate hbv-artic    # change to activate the custom environment

# Run the pipeline
nextflow run main.nf -profile hbv,slurm -resume
