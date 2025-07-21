#!/bin/bash
#SBATCH -J CR_BM2  #job name
#SBATCH --partition=amd        #Partition/queue to run job #GPU
#SBATCH -t 6:00:00               #max time 3 hours
#SBATCH --cpus-per-task=32        #nr os cpu requested
#SBATCH --mem=256GB 

CONFIG_CSV=/gpfs/helios/home/rostamne/CVS_SplVacc_BM__Ghr_Thy_Jul_24_scRNA/code/BM2/config.csv

export PATH=/gpfs/helios/home/rostamne/msc/software/cellranger-9.0.1/:$PATH

cellranger multi \
  --id=GEX_BM1_multi \
  --csv=$CONFIG_CSV \
  --localcores=32 \
  --localmem=256





