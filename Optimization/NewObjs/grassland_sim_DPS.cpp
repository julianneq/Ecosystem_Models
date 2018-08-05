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
#define reliability_threshold 0.85

double h, I, F_x, S, F_y;

double phi = 0.2;
double sigma = 5.0;

double r_x = 4000.0;
double c2 = 3.15;
double G0 = 500.0;
double d3 = 600.0;
double c3 = d3 - G0;
double i0 = 0.1;
double alpha = 0.122/0.21;
double L_x = 1.0;
double beta = 0.01/0.21;
double d1 = 500.0;
double c1 = d1 + i0*d1/(1-i0);
double x_in = 1.0;

double r_y = 3000.0;
double r_S = 2000.0;
double d4 = 500.0;
double K = 2.0;
double p = 0.1;
double L_y = 1.0;
double gamma_param = 1.0/alpha;
double delta = beta/alpha;

double dt = 1.0;

double x_0 = 2000.0;
double y_0 = 2000.0;

double u = 0.0;

int nvars = 3;
int nobjs = 2;
int nconstrs = 0;

namespace ublas = boost::numeric::ublas;
namespace tools = boost::math::tools;
using namespace std;

ublas::vector<double> min_annual_x(N);

ublas::vector<double> x(T+3); // grass biomass
ublas::vector<double> y(T+3); // woody biomass
ublas::vector<double> dxdt(T+2);
ublas::vector<double> dydt(T+2);
ublas::vector<double> Q(T+2);
ublas::vector<double> z(T+2);
ublas::vector<double> a(T);

void grassland_problem(double* vars, double* objs, double* constrs) 
{
  // find stocking rate
  h = (int) vars[0] + 0.5;

  // initialize variables
  zero(min_annual_x);
  zero(x);
  zero(y);
  zero(dxdt);
  zero(dydt);
  zero(Q);
  zero(z);
  zero(a);

  // run lake model simulation
  for (int i = 0; i < N; i++){
    ublas::vector<double> eps = generateSOW(T+3, 0.0, sigma);

    x[0] = x_0 + eps[0];
    y[0] = y_0 - eps[0];

    //implement the grassland model from Carpenter et al. 2015
    for (int t = 0; t < T+2; t++){
      I = (x[t] + c1*i0)/(x[t] + c1);
      F_x = 1/(x[t] + alpha*y[t] + beta);
      Q[t] = max(0.0, c2*(1 - G0/x[t])/(x[t] + c3));
      S = r_S*(1 - x[t]/(K*(x[t] + d4)));
      F_y = 1/(y[t] + gamma_param*x[t] + delta);
      dxdt[t] = r_x*I*F_x*x[t] - L_x*x[t] + x_in;
      dydt[t] = r_y*I*F_y*y[t] + S*y[t]/(y[t] + p);
      if (t < 2){ // just add noise
        z[t] = h*x[t]*Q[t];
        x[t+1] = x_0 + eps[t+1];
        y[t+1] = y_0 - eps[t+1];
      }
      else{
        z[t] = h*Q[t]*x[t] - phi*h*Q[t-1]*x[t-1];
        x[t+1] = max(1.0, x[t] + phi*(x[t] - x[t-1]) + dxdt[t]*dt - phi*dxdt[t-1]*dt
          - max(0.1, z[t])*dt + eps[t+1]);
        a[t-2] = L_y + min(1.0, max(0.0, pow(abs(x[t]-vars[1])/vars[2],3.0)));
        y[t+1] = max(1.0, y[t] + dydt[t]*dt - a[t-2]*y[t]*dt + u*a[t-2]*y[t-1]);
        if (x[t+1] < min_annual_x[i]){
          min_annual_x[i] = x[t+1];
        }
      }
    }
  }

  objs[0] = -h; // stocking rate
  objs[1] = -1*vsum(min_annual_x)/N; // average minimum grassland biomass

  //constrs[0] = max(0.0, reliability_threshold - (-1*objs[1]));
  
  min_annual_x.clear();
  x.clear();
  y.clear();
  dxdt.clear();
  dydt.clear();
  z.clear();
  a.clear();
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
  BORG_Problem problem = BORG_Problem_create(nvars, nobjs, nconstrs, grassland_problem);

  // Set all the parameter bounds and epsilons
  BORG_Problem_set_bounds(problem, 0, 1.0, 1000.0);
  BORG_Problem_set_bounds(problem, 1, 0.0, 3500.0);
  BORG_Problem_set_bounds(problem, 2, 0.0, 3500.0);

  BORG_Problem_set_epsilon(problem, 0, 1.0); // stocking rate
  BORG_Problem_set_epsilon(problem, 1, 1.0); // average grassland biomass

  //This is set up to run only one seed at a time
  char outputFilename[256];
  char runtime[256];
  FILE* outputFile = NULL;
  sprintf(outputFilename, "./sets/GrasslandDPS_S%d.set", seed);
  sprintf(runtime, "./runtime/GrasslandDPS_S%d.runtime", seed);

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
