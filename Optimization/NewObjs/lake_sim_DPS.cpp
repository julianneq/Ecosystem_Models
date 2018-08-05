/* main-lake_DPS.cpp
   
  Riddhi Singh, May, 2014
  The Pennsylvania State University
  rus197@psu.edu

  Adapted by Tori Ward, July 2014 
  Cornell University
  vlw27@cornell.edu

  Adapted by Jonathan Herman and David Hadka, Sept-Dec 2014
  Cornell University and The Pennsylvania State University

  Adapted by Julianne Quinn, July 2015 as DPS problem
  Cornell University
  jdq8@cornell.edu

  A multi-objective represention of the lake model from Carpenter et al., 1999
  This simulation is designed for optimization with either Borg or the MOEAFramework.

  Stochasticity is introduced by natural phosphorous inflows. 
  These follow a lognormal distribution with specified mu and sigma.

  Decision variable
    vars : vector of length 3n describing the policy parameters
    n is the number of RBFs, each with a center, radius and weight
    format of vars is [c1, r1, w1, c2, r2, w2, ..., cn, rn, wn]
    all weights must be between 0 and 1 and sum to 1
    The weights determined by Borg aren't constrained to sum to 1,
    but are normalized to do so in the model simulation
    The actual decision, the amount of pollution emitted by the town, is Y,
    a function of vars and the current state of the lake (lake_state)

  Objectives
  1: max avg P concentration
  2: mean economic benefits
  3: mean inertia
  4: mean reliability

  Constraints
  Reliability must be > 85%

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sstream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/tools/roots.hpp>
#include "borg/moeaframework.h"
#include "boostutil.h"
#include "borg/borgms.h"
#include "generate-scenarios.h"

#define T 100
#define N 100

double phi = 0.1;
double sigma = 0.35;

double alpha = 0.4;
double delta = 0.98;

double h = 1 - exp(-0.29);
double s = 1 - h;
double r = 0.019;
double m = 4.0;
double q = 4.0;
double b = 0.002;

double dt = 1.0;

double x_0 = 1.0;
double y_0 = 330.0;

double u = 0.0;

int nvars = 4;
int nobjs = 2;
int nconstrs = 0;

namespace ublas = boost::numeric::ublas;
namespace tools = boost::math::tools;
using namespace std;

ublas::vector<double> discounted_benefit(N);
ublas::vector<double> max_annual_P(N);

ublas::vector<double> x(T+3);
ublas::vector<double> y(T+3);
ublas::vector<double> dxdt(T+2);
ublas::vector<double> dydt(T+2);
ublas::vector<double> a(T);

void lake_problem(double* vars, double* objs, double* constrs) 
{
  // initialize variables
  zero(discounted_benefit);
  zero(max_annual_P);
  zero(x);
  zero(y);
  zero(dxdt);
  zero(dydt);
  zero(a);

  // run lake model simulation
  for (int i = 0; i < N; i++){
    ublas::vector<double> eps = generateSOW(T+3, 0.0, sigma);

    x[0] = x_0 + eps[0];
    y[0] = y_0 - eps[0];

    //implement the lake model from Carpenter et al. 2015
    for (int t = 0; t < T+2; t++){
      dxdt[t] = -(s + h)*x[t] + r*y[t]*pow(x[t],q) / (pow(m,q) + pow(x[t],q));
      dydt[t] = s*x[t] - b*y[t] - r*y[t]*pow(x[t],q) / (pow(m,q) + pow(x[t],q));
      if (t < 2){
        x[t+1] = x_0 + eps[t+1];
        y[t+1] = y_0 - eps[t+1];
      }
      else{
        a[t-2] = min(2.5, max(max(vars[0]*x[t] + vars[1], vars[2]*x[t] + vars[3]), 0.1));
        x[t+1] = max(0.1, x[t] + (u + phi)*(x[t] - x[t-1]) - u*phi*(x[t-1] - x[t-2])
          + dxdt[t]*dt - (u + phi)*dxdt[t-1]*dt + u*phi*dxdt[t-2]*dt
          + (1 - (u + phi) + u*phi)*a[t-2] + eps[t+1]);
        y[t+1] = max(1.0, y[t] + dydt[t]*dt);
        discounted_benefit[i] += alpha*a[t-2]*pow(delta,t-2);
        if (x[t+1] > max_annual_P[i]){
        	max_annual_P[i] = x[t+1];
        }
      }
    }
  }

  objs[0] = -1*vsum(discounted_benefit)/N; // average economic benefit
  objs[1] = -1*vsum(max_annual_P)/N; // average maximum P concentration

  //constrs[0] = max(0.0, reliability_threshold - (-1*objs[1]));
  
  discounted_benefit.clear();
  max_annual_P.clear();
  x.clear();
  y.clear();
  dxdt.clear();
  dydt.clear();
  a.clear();

}

double root_function(double x) {
  return pow(x,q)/(1+pow(x,q)) - b*x;
}

bool root_termination(double min, double max) {
  return abs(max - min) <= 0.000001;
}

int main(int argc, char* argv[]) 
{
  // setting random seed
  unsigned int seed = atoi(argv[1]);
  srand(seed);
  int NFE = atoi(argv[2]);

  // interface with Borg-MS
  BORG_Algorithm_ms_startup(&argc, &argv);
  BORG_Algorithm_ms_max_evaluations(NFE);
  BORG_Algorithm_output_frequency(NFE/1000);

  // Define the problem with decisions, objectives, constraints and the evaluation function
  BORG_Problem problem = BORG_Problem_create(nvars, nobjs, nconstrs, lake_problem);

  // Set all the parameter bounds and epsilons
  BORG_Problem_set_bounds(problem, 0, -5.0, -0.2);
  BORG_Problem_set_bounds(problem, 1, 0.0, 10.0);
  BORG_Problem_set_bounds(problem, 2, 0.2, 5.0);
  BORG_Problem_set_bounds(problem, 3, -10.0, 0.0);

  BORG_Problem_set_epsilon(problem, 0, 0.01); // average economic benefit
  BORG_Problem_set_epsilon(problem, 1, 0.0001); // average reliability

  //This is set up to run only one seed at a time
  char outputFilename[256];
  char runtime[256];
  FILE* outputFile = NULL;
  sprintf(outputFilename, "./sets/LakeDPS_S%d.set", seed);
  sprintf(runtime, "./runtime/LakeDPS_S%d.runtime", seed);

  BORG_Algorithm_output_runtime(runtime);

  BORG_Random_seed(seed);
  BORG_Archive result = BORG_Algorithm_ms_run(problem); // this actually runs the optimization

  //If this is the master node, print out the final archive
  if (result != NULL){
    outputFile = fopen(outputFilename, "w");
    if(!outputFile){
      BORG_Debug("Unable to open final output file\n");
    }
    BORG_Archive_print(result, outputFile);
    BORG_Archive_destroy(result);
    fclose(outputFile);
  }

  BORG_Algorithm_ms_shutdown();
  BORG_Problem_destroy(problem);

  return EXIT_SUCCESS;

}
