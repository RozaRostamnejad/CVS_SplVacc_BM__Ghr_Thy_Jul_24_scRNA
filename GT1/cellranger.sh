#!/bin/bash
#SBATCH -J CR_GT1_new
#SBATCH --partition=amd
#SBATCH -t 10:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=256GB

CONFIG_CSV=/gpfs/helios/home/rostamne/CVS_SplVacc_BM__Ghr_Thy_Jul_24_scRNA/code/GT1/config.csv


cellranger multi \
  --id=GEX_GT1_new_multi \
  --csv=$CONFIG_CSV \
  --localcores=32 \
  --localmem=256
