#!/bin/bash
#SBATCH -J CR_BM2
#SBATCH --partition=amd
#SBATCH -t 10:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=256GB

cellranger multi \
  --id=GEX_BM2_new_multi \
  --csv=/gpfs/helios/home/rostamne/CVS_SplVacc_BM__Ghr_Thy_Jul_24_scRNA/code/BM2/config.csv \
  --localcores=32 \
  --localmem=256
