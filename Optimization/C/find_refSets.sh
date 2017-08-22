#/bin/bash

source /etc/profile.d/modules.sh
module load python-2.7.5

python pareto.py ./sets/*.set -o 4-5 -e 0.01 0.001 --output Lake_DPS.resultfile --delimiter=" " --comment="#"
cut -d ' ' -f 5-6 Lake_DPS.resultfile >Lake_DPS.reference