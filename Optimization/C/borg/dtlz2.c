/* Copyright 2012-2013 The Pennsylvania State University
 *
 * This software was written by David Hadka and others.
 * 
 * The use, modification and distribution of this software is governed by the
 * The Pennsylvania State University Research and Educational Use License.
 * You should have received a copy of this license along with this program.
 * If not, contact <preed@engr.psu.edu>.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "moeaframework.h"

#define PI 3.14159265358979323846

int nvars;
int nobjs;

void evaluate(double* vars, double* objs) {
	int i;
	int j;
	int k = nvars - nobjs + 1;
	double g = 0.0;

	for (i=nvars-k; i<nvars; i++) {
		g += pow(vars[i] - 0.5, 2.0);
	}

	for (i=0; i<nobjs; i++) {
		objs[i] = 1.0 + g;

		for (j=0; j<nobjs-i-1; j++) {
			objs[i] *= cos(0.5*PI*vars[j]);
		}

		if (i != 0) {
			objs[i] *= sin(0.5*PI*vars[nobjs-i-1]);
		}
	}
}

int main(int argc, char* argv[]) {
	if (argc <= 1) {
		nobjs = 2;
		nvars = 11;
	} else if (argc == 2) {
		nobjs = atoi(argv[1]);
		nvars = nobjs + 9;
	} else if (argc >= 3) {
		nobjs = atoi(argv[2]);
		nvars = atoi(argv[1]);
	} 

	double vars[nvars];
	double objs[nobjs];

#ifdef USE_SOCKET
	MOEA_Init_socket(nobjs, 0, NULL);
#else
	MOEA_Init(nobjs, 0);
#endif

	while (MOEA_Next_solution() == MOEA_SUCCESS) {
		MOEA_Read_doubles(nvars, vars);
		evaluate(vars, objs);
		MOEA_Write(objs, NULL);
	}

	return EXIT_SUCCESS;
}

