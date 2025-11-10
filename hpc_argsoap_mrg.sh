#!/bin/bash

#SBATCH --job-name="AMR_HMR"
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=40G
#SBATCH --array=1-18%4
#SBATCH --output=/ibiscostorage/dcorso/ideARG_2020/heavy_metal_res_READBASED/logs_slurm/slurm-%A_%a.out
#SBATCH --partition=parallel

#
# Copyright 2023 Davide Corso
#

# CREATE THE ENVIRONMENT
# mamba create -n args_oap -c bioconda -c conda-forge args_oap


#### PRIMA DI ESEGUIRE, ATTIVARE AMBIENTE CONDA E CREARE DB con: args_oap make_db -i <metal_sequences.fasta>


singleSample="$(tail -n +$SLURM_ARRAY_TASK_ID /ibiscostorage/dcorso/ideARG_2020/idearg_samples.txt | head -n1)"

file_database="/ibiscostorage/dcorso/ideARG_2020/heavy_metal_res_READBASED/args_oap_METALDB/metal_sequences.fasta"
file_structure="/ibiscostorage/dcorso/ideARG_2020/heavy_metal_res_READBASED/args_oap_METALDB/structure.txt"

echo "SAMPLE:  ----  $singleSample"


folder_reads="/ibiscostorage/dcorso/ideARG_2020/trimmed_reads_from_kbase/deinterleave/${singleSample}"
folder_output="/ibiscostorage/dcorso/ideARG_2020/heavy_metal_res_READBASED/output/${singleSample}/"

mkdir -p "${folder_output}"

args_oap stage_one -i $folder_reads -o $folder_output -f fq -t 50 --database $file_database
args_oap stage_two -i $folder_output -t 50 --database $file_database --structure1 $file_structure
