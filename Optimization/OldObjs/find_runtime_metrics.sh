NSEEDS=50
SEEDS=$(seq 1 ${NSEEDS})
JAVA_ARGS="-cp MOEAFramework-2.4-Demo.jar"

for SEED in ${SEEDS}
do
	NAME=Runtime_Metrics_S${SEED}
	PBS="\
	#PBS -N ${NAME}\n\
	#PBS -l nodes=1\n\
	#PBS -l walltime=1:00:00\n\
	#PBS -o output/${NAME}\n\
	#PBS -e error/${NAME}\n\
	cd \$PBS_O_WORKDIR\n\
	java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.ResultFileEvaluator \
		-d 2 -i ./objs/GrasslandDPS_S${SEED}.obj -r Grassland_DPS.reference \
		-o ./metrics/GrasslandDPS_S${SEED}.metrics"
	echo -e $PBS | qsub
done

# for SEED in ${SEEDS}
# do
# 	NAME=Runtime_Metrics_S${SEED}
# 	PBS="\
# 	#PBS -N ${NAME}\n\
# 	#PBS -l nodes=1\n\
# 	#PBS -l walltime=1:00:00\n\
# 	#PBS -o output/${NAME}\n\
# 	#PBS -e error/${NAME}\n\
# 	cd \$PBS_O_WORKDIR\n\
# 	java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.ResultFileEvaluator \
# 		-d 2 -i ./objs/LakeDPS_S${SEED}.obj -r Lake_DPS.reference \
# 		-o ./metrics/LakeDPS_S${SEED}.metrics"
# 	echo -e $PBS | qsub
# done
