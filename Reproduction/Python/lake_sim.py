from __future__ import division
import numpy as np
from computeSpectra import computeSpectra
from make_plots import plotTimeSeries, plotSpectra, plotCrash

# Model Parameters
phi = 0.1 # auto-regressive coefficient of random natural P inputs
sigma = 0.35 # standard deviation of random natural P inputs

h = 1 - np.exp(-0.29) # P hydraulic washout rate
s = 1 - h # P sedimentation rate
r = 0.019 # max P recycling rate
m = 4 # half-saturation coefficient for P recycling
q = 4 # shape parameter of P recycling curve
b = 0.002 # P loss rate from sediment

# Simulation parameters
dt = 1 # length of time step in years

# Initial conditions
x0 = 1 # initial lake P concentration
y0 = 330 # initial sediment P concentration

# Control parameter determining anthropogenic inputs, U(B)
# U(B)=1-u*B where U(B) is the manager's intervention and B is the backshift operator
# u = 0.6 corresponds to a manager trying to minimize short-term variance
# u = 0.0 corresponds to a manager ignoring variance in their management strategy
u = np.array([0.6, 0.0])

def lake_sim(llambda, T, eps, burnin):
    # State variables
    x = np.zeros([len(u),T])
    y = np.zeros([len(u),T])
    dxdt = np.zeros([T-1])
    dydt = np.zeros([T-1])
    
    # anthropogenic emission under noise-management policy
    a = np.zeros([T-3])
    
    # simulate remaining time steps for both management strategies
    # set random perturbation of initial conditions
    x[:,0] = x0 + eps[0]
    y[:,0] = y0 - eps[0]
    for i in range(len(u)):
        for t in range(T-1):
            dxdt[t] = -(s+h)*x[i,t] + (r*y[i,t]*x[i,t]**q / (m**q + x[i,t]**q))
            dydt[t] = s*x[i,t] - b*y[i,t] - (r*y[i,t]*x[i,t]**q / (m**q + x[i,t]**q))
            if t < 2: # just add random natural P inputs for first 2 time steps (transferred from sediment to lake)
                x[i,t+1] = x0 + eps[t+1]
                y[i,t+1] = y0 - eps[t+1]
            else: # add anthropogenic + random natural P inputs + P recycled from sediment
                x[i,t+1] = max(0.1, x[i,t] + (u[i]+phi)*(x[i,t]-x[i,t-1]) - u[i]*phi*(x[i,t-1]-x[i,t-2]) \
                    + dxdt[t]*dt - (u[i] + phi)*dxdt[t-1]*dt + u[i]*phi*dxdt[t-2]*dt \
                    + (1-(u[i]+phi)+u[i]*phi)*llambda + eps[t+1])
                y[i,t+1] = max(1.0, y[i,t] + dydt[t]*dt)
                if u[i] != 0.0: # mean emission + difference between noise-management and no noise-mangement policies
                    a[t-2] = llambda + x[i,t+1] - max(0.1, x[i,t] + phi*(x[i,t]-x[i,t-1]) \
                    + dxdt[t]*dt - phi*dxdt[t-1]*dt + (1-phi)*llambda + eps[t+1])
                
    # remove burn-in
    x = x[:,burnin::]
    y = y[:,burnin::]
    a = a[(burnin-3)::]
    
    return (x, y, a)
    
# run simulation
llambda = 0.5 # mean rate of anthropogenic P inputs
burnin = 500
T = 2**11 + burnin # total number of time steps
#eps = np.random.normal(0,sigma,T) # generate noise
eps = np.loadtxt('./../R/epsilon_lake.txt') # read in to verify same output from R and Python
x, y, a = lake_sim(llambda, T, eps, burnin)

# plot time series of x w/ and w/o variance control and compare with R
labels = ['Var Mgmt - R','Var Mgmt - Python','Nominal - R','Nominal - Python']
colors = ['#a50f15','#fb6a4a','#08519c','#6baed6']  
plotTimeSeries(x, len(u), T, burnin, labels, colors, './../R/lakeP.txt','Lake/LakeTimeSeries')

# compute spectra
k = 30 # number of Slepian windows
NW = k/2 # time half bandwith parameter
Slep_real_sq, weights = computeSpectra(x, NW, k, burnin, T)

# plot spectra and compare with R
plotSpectra(len(u), T, burnin, Slep_real_sq, weights, labels, colors, './../R/lake_spectra.txt','Lake/LakeSpectra')

# run simulation for multiple llambdas
llambdas = np.arange(0.25, 2.5+2.25/14.0, 2.25/14.0)
T = 5500 # total number of time steps
#eps = np.random.normal(0,sigma,T) # generate noise
eps = np.loadtxt('./../R/epsilon2_lake.txt') # read in to verify same output from R and Python
x = np.zeros([len(llambdas),len(u),T-burnin])
y = np.zeros([len(llambdas),len(u),T-burnin])
a = np.zeros([len(llambdas), T-burnin])
for i in range(len(llambdas)):
    x[i,:,:], y[i,:,:], a[i,:] = lake_sim(llambdas[i], T, eps, burnin)

# find fraction of days in oligotrophic state
pct_oligo = (x <= 3.0).sum(axis=2)/(T-burnin)

# plot pct_oligo and compare with R
plotCrash(len(u), llambdas, pct_oligo, labels, colors, './../R/pct_oligo.txt','Lake/PctOligo')