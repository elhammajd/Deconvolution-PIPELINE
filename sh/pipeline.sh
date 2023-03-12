#!/bin/bash
#SBATCH --account=def-ubcxzh
#SBATCH --output=./output/%x-%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=elhmajd@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH --ntasks=2  
#SBATCH --mem-per-cpu=16G # 2GiB of memery 
#SBATCH -t 00-01:20
Rscript PIPELINE.R
