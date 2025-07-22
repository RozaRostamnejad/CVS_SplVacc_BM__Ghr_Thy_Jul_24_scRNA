#!/bin/bash
#SBATCH -J CR_GT2  #job name
#SBATCH --partition=amd        #Partition/queue to run job #GPU
#SBATCH -t 6:00:00               #max time 3 hours
#SBATCH --cpus-per-task=32        #nr os cpu requested
#SBATCH --mem=256GB 

CONFIG_CSV=/gpfs/helios/home/rostamne/CVS_SplVacc_BM__Ghr_Thy_Jul_24_scRNA/code/GT2/config.csv


cellranger multi \
  --id=GEX_GT2_multi \
  --csv=$CONFIG_CSV \
  --localcores=32 \
  --localmem=256





