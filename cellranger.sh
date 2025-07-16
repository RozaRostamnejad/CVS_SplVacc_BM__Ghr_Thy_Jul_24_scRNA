#!/bin/bash
#SBATCH -J cellranger_BM1  #job name
#SBATCH --partition=amd        #Partition/queue to run job #GPU
#SBATCH -t 3:00:00               #max time 3 hours
#SBATCH --cpus-per-task=32        #nr os cpu requested
#SBATCH --mem=156GB 

CONFIG_CSV=/gpfs/helios/home/rostamne/CVS_SplVacc_BM__Ghr_Thy_Jul_24_scRNA/code/config.csv

export PATH=/gpfs/helios/home/rostamne/msc/software/cellranger-9.0.1/:$PATH

cellranger multi \
  --id=GEX_BM1_multi \
  --csv=$CONFIG_CSV \
  --localcores=32 \
  --localmem=156





