#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=16G
#SBATCH --time=0-00:30:00
#SBATCH --output=log.dat
#SBATCH --mail-user=ahmedi@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="DarkMatterHalos Trial 1"

# Put commands for executing job below this line
# This example is loading Python 2.7.8 and then
# writing out the version of Python
setpkgs -a python2.7.8_intel14
python --version
python multicore_example.py
