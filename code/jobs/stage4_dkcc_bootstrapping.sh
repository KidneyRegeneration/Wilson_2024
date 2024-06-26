#!/bin/bash

### slurm job options: line must begin with '#SBATCH'

#SBATCH --job-name=any_name    # job name
#SBATCH --mail-type=END,FAIL    # mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=your.name@sund.ku.dk    # email address to receive the notification    
#SBATCH -c 2    # number of requested cores
#SBATCH --mem=8gb    # total requested RAM
##SBATCH --gres=gpu:1    # request 1 gpu; if you need it, copy this line and delete one '#', otherwise please do not include this line
#SBATCH --time=0-00:05:00               # max. running time of the job, format in D-HH:MM:SS
#SBATCH --output=any_name_%j.log   # Standard output and error log, '%j' gives the job ID

### write your own scripts below 

module load dangpu_libs python/3.7.13     
echo "I am running a job with the slurm job script"
python /this/is/a/fake/path/empty_python_script.py