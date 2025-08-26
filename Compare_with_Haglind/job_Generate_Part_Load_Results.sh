#!/bin/bash
#SBATCH -p long # partition (queue)
#SBATCH -N 1 # (leave at 1 unless using multi-node specific code)
#SBATCH -n 1 # number of cores
#SBATCH --mem-per-cpu=160000 # memory per core
#SBATCH --job-name="Run_Part_load"
#SBATCH --mail-user=ib011@bucknell.edu # address to email
#SBATCH --mail-type=ALL # mail events (NONE, BEGIN, END, FAIL, ALL)
module load matlab
matlab -singleCompThread -nodisplay -nosplash -r Generate_Part_Load_Results



