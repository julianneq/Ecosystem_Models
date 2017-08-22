#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1
#PBS -j oe

source /etc/profile.d/modules.sh
module load python-2.7.5
cd /home/fs02/pmr82_0001/jdq8/Ecosystem_Models/Optimization
python lake_sim.py