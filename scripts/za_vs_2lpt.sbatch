#!/bin/bash
#SBATCH --nodes=6
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-01:00:00
#SBATCH --output=results/za_vs_2lpt.dat
#SBATCH --mail-user=ahmedi@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="DarkMatterHalos Trial 3"

# Put commands for executing job below this line
# This example is loading Python 2.7.8 and then
# writing out the version of Python
setpkgs -a python2.7.8
python --version
python za_vs_2lpt.py 62 box1
python za_vs_2lpt.py 62 box2
python za_vs_2lpt.py 62 box3
