import numpy as np
import os
from sys import *
from borg import *
from math import *

# Decision variable
# vars : vector of length 4 describing the policy parameters
# format of vars is [m1, b1, m2, b2]
# where m1 and b1 are the slope and intercept of the descending line
# and m2 and b2 are the slope and intercept of the ascending line

# Model Parameters
phi = 0.1 # auto-regressive coefficient of random natural P inputs
sigma = 0.35 # standard deviation of random natural P inputs

alpha = 0.4
delta = 0.98 # discount rate

h = 1 - np.exp(-0.29) # P hydraulic washout rate
s = 1 - h # P sedimentation rate
r = 0.019 # max P recycling rate
m = 4.0 # half-saturation coefficient for P recycling
q = 4.0 # shape parameter of P recycling curve
b = 0.002 # P loss rate from sediment

# Simulation parameters
dt = 1.0 # length of time step in years

# Initial conditions
x0 = 1.0 # initial lake P concentration
y0 = 330.0 # initial sediment P concentration

# Control parameter determining anthropogenic inputs, U(B)
# U(B)=1-u*B where U(B) is the manager's intervention and B is the backshift operator
# u = 0.6 corresponds to a manager trying to minimize short-term variance
# u = 0.0 corresponds to a manager ignoring variance in their management strategy
u = 0.0

# number of time steps (T) and ensemble members (S)
T = 100
S = 100

# set optimization parameters
nvars = 4 # number of decision variables (parameters of RBF policy)
nobjs = 2 # number of objectives (1: average economic benefits, 2: average reliability)
nconstrs = 1 # average reliability >= reliability_threshold
reliability_threshold = 0.85

def lake_sim(*vars):
    # Initialize arrays to store objective function values and constraints
    objs = [0.0]*nobjs
    constrs = [0.0]*nconstrs
    
    #Initialize arrays to store performance of policy in each realization
    discounted_benefit = np.zeros([S])
    yrs_pcrit_met = np.zeros([S])
    
    # generate S ensemble of T+3 years of Gaussian random noise
    eps = np.zeros([S,T+3])
    for i in range(S):
        eps[i,:] = np.random.normal(0,sigma,T+3)
    
    # run lake model simulation
    for i in range(S):
        # State variables
        x = np.zeros([T+3])
        y = np.zeros([T+3])
        dxdt = np.zeros([T+2])
        dydt = np.zeros([T+2])
        
        # anthropogenic P emissions
        a = np.zeros([T])
        
        # run simulatiobn
        # set random perturbation of initial conditions for fist time step
        x[0] = x0 + eps[i,0]
        y[0] = y0 - eps[i,0]
        for t in range(T+2):
            dxdt[t] = -(s+h)*x[t] + (r*y[t]*x[t]**q / (m**q + x[t]**q))
            dydt[t] = s*x[t] - b*y[t] - (r*y[t]*x[t]**q / (m**q + x[t]**q))
            if t < 2: # just add random natural P inputs for first 2 time steps
                x[t+1] = x0 + eps[i,t+1]
                y[t+1] = y0 - eps[i,t+1]
            else: # add anthropogenic + random natural P inputs + P recycled from sediment
                a[t-2] = min(1.5, max(max(vars[0]*x[t] + vars[1], vars[2]*x[t] + vars[3]), 0.01))
                x[t+1] = max(0.1, x[t] + (u + phi)*(x[t] - x[t-1]) - u*phi*(x[t-1] - x[t-2]) \
                    + dxdt[t]*dt - (u + phi)*dxdt[t-1]*dt + u*phi*dxdt[t-2]*dt \
                    + (1 - (u + phi) + u*phi)*a[t-2] + eps[i,t+1])
                y[t+1] = max(1.0, y[t] + dydt[t]*dt)
                discounted_benefit[i] += alpha*a[t-2]*delta**(t-2)
                
                if x[t+1] <= 3.0: # critical P threshold
                    yrs_pcrit_met[i] += 1
            
    # compute objectives and constraints
    objs[0] = -1*np.sum(discounted_benefit)/S
    objs[1] = -1*np.sum(yrs_pcrit_met)/(T*S)
    constrs[0] = max(0.0, reliability_threshold - (-1*objs[1]))
    
    return (objs, constrs)

# Set up interface with Borg for optimization
borg = Borg(nvars, nobjs, nconstrs, lake_sim)
borg.setBounds(*[[-5.0,-0.2], [0.0,10.0], [0.2,5.0], [-10.0,0.0]])
borg.setEpsilons(0.01, 0.0001)

result = borg.solve({"maxEvaluations":1000})

f = open(os.getcwd() + '/sets/Lake/Borg_Lake_DPS_S' + \
    str(j+1) + '.set','w')

f.write('# Borg Optimization Results\n')
f.write('# First ' + str(nvars) + ' are the decision variables, ' \
    'last ' + str(nobjs) + ' are the objective values\n')

for solution in result:
    line = ''
    for i in range(len(solution.getVariables())):
        line = line + (str(solution.getVariables()[i])) + ' '

    for i in range(len(solution.getObjectives())):
        line = line + (str(solution.getObjectives()[i])) + ' '

    f.write(line[0:-1]+'\n')

f.write("#")
    
f.close()
