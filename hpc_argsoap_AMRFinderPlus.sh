#!/bin/bash

#SBATCH --job-name="AMR_Soap"
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=40G
#SBATCH --array=1-18%5
#SBATCH --output=/ibiscostorage/dcorso/ideARG_2020/args_oap_pipeline/logs_slurm/slurm-%A_%a.out
#SBATCH --partition=parallel

#
# Copyright 2023 Davide Corso
#

# AMRFinderPlus on Reads using args_oap tool

# CREATE THE ENVIRONMENT
# mamba create -n args_oap -c bioconda -c conda-forge args_oap


singleSample="$(tail -n +$SLURM_ARRAY_TASK_ID /ibiscostorage/dcorso/ideARG_2020/idearg_samples.txt | head -n1)"

amrfp_database="/ibiscostorage/dcorso/ideARG_2020/args_oap_pipeline/argsoap_db/sequences.fasta"
amrfp_structure="/ibiscostorage/dcorso/ideARG_2020/args_oap_pipeline/argsoap_db/structure.txt"

echo "SAMPLE:  ----  $singleSample"


folder_reads="/ibiscostorage/dcorso/ideARG_2020/trimmed_reads_from_kbase/deinterleave/${singleSample}"
folder_output="/ibiscostorage/dcorso/ideARG_2020/args_oap_pipeline/output/${singleSample}/"

mkdir -p "${folder_output}"

args_oap stage_one -i $folder_reads -o $folder_output -f fq -t 50 --database $amrfp_database
args_oap stage_two -i $folder_output -t 50 --database $amrfp_database --structure1 $amrfp_structure
