#/bin/bash

source /etc/profile.d/modules.sh
module load python-2.7.5

# python pareto.py ./sets/Lake*.set -o 4-5 -e 0.01 0.001 --output Lake_DPS.resultfile --delimiter=" " --comment="#"
# cut -d ' ' -f 5-6 Lake_DPS.resultfile >Lake_DPS.reference

python pareto.py ./sets/Grassland*.set -o 3-4 -e 1.0 1.0 --output Grassland_DPS.resultfile --delimiter=" " --comment="#"
cut -d ' ' -f 4-5 Grassland_DPS.resultfile >Grassland_DPS.reference