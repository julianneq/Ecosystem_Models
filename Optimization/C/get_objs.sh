#!/bin/bash
NSEEDS=50
SEEDS=$(seq 1 ${NSEEDS})

for SEED in ${SEEDS}
do
	awk 'BEGIN {FS=" "}; /^#/ {print $0}; /^[^#/]/ {printf("%s %s\n",$5,$6)}' ./runtime/LakeDPS_S${SEED}.runtime \
	 	>./objs/LakeDPS_S${SEED}.obj
done