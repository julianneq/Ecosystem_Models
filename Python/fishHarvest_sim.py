from __future__ import division
import numpy as np
from computeSpectra import computeSpectra
from make_plots import plotTimeSeries, plotSpectra, plotCrash

# Model parameters
phi = 0.4 # auto-regressive coefficient of noise
sigma = 0.1 # standard deviation of noise

r = 0.75 # logistic growth parameter
c = 1.35 # sigmoid predation consumption coefficient
b = 1 # half-saturation constant
q = 2 # shape parameter for logistic growth
K = 8.8  # critical K logistic growth parameter is about 6.5 for h=0, 7.6 for h=0.05, 9 for h=0.1
h = 0.02 # harvest coefficient: fraction of fish population harvested

# Initial conditions
x0 = (r-h)*(K/r) # initial fish population: deterministic equilibrium for logistic growth

def fishHarvest_sim(phi, h, u, T, eps, burnin, control, crashCount):
    # State variables
    x = np.zeros([len(u),T]) # fish population
    z = np.zeros([len(u),T]) # fish harvest
    dx = np.zeros([T])
    g = np.zeros([T-1])
    
    # simulate remaining time steps for both management strategies
    # set random perturbation of initial conditions
    x[:,0] = x0 + eps[0]
    z[:,0] = h*x[:,0]
    for i in range(len(u)):
        for t in range(T-1):
            g[t] = r*x[i,t]*(1-x[i,t]/K) - c*x[i,t]**q / (b**q + x[i,t]**q) # logistic growth
            if t < 2: # just add random noise
                x[i,t+1] = x0 + eps[t+1]
                z[i,t+1] = h*x[i,t+1] # harvest for next time step
            else:
                if control == 'harvest':
                    x[i,t+1] = max(0.1, x[i,t] + phi*dx[t] + g[t] - phi*g[t-1] - z[i,t] + eps[t+1])
                    z[i,t+1] = h*x[i,t+1] - (u[i]+phi)*h*x[i,t] + u[i]*phi*h*x[i,t-1]
                else: # control = 'noise'
                    x[i,t+1] = max(0.1, x[i,t] + (u[i]+phi)*dx[t] - u[i]*phi*dx[t-1] + \
                        (g[t]-z[i,t]) - (u[i]+phi)*(g[t-1]-z[i,t-1]) + u[i]*phi*(g[t-2]-z[i,t-2]) + eps[t+1])
                    z[i,t+1] = h*x[i,t+1]
                    
                if crashCount == True and x[i,t+1] <= 0.12: # re-set population if finding average time to crash
                    x[i,t+1] = x0
            
            dx[t+1] = x[i,t+1] - x[i,t]
    
    # remove burn-in
    x = x[:,burnin::]
    z = z[:,burnin::]
    
    return (x, z)
    
def runSimulation(u, n, labels, colors, control):
    # define inputs
    burnin = 100
    T = 2**11 + burnin # total number of time steps
    #eps = np.random.normal(0,sigma,T) # generate noise
    eps = np.loadtxt('./../R/epsilon_fish.txt') # read in to verify same output from R and Python
    
    # run simulation
    x, z = fishHarvest_sim(phi, h, u, T, eps, burnin, control, False)
    if control == 'noise': # run white noise simulation (phi = 0.0 and u = 0.0)
        x2, z2 = fishHarvest_sim(0.0, h, [0.0], T, eps, burnin, control, False)
        z = np.concatenate((z,z2),0)
    
    # plot time series of z and compare with R
    plotTimeSeries(z, n, T, burnin, labels, colors, './../R/fishHarvest_' + control + '_mgmt.txt', \
        'Fish/HarvestTimeSeries_' + control + '_mgmt')
        
    # compute spectra
    k = 40 # number of Slepian windows
    NW = k/2 # time half bandwith parameter
    Slep_real_sq, weights = computeSpectra(z, NW, k, burnin, T)
    
    # plot spectra and compare with R
    plotSpectra(n, T, burnin, Slep_real_sq, weights, labels, colors, \
        './../R/fishSpectra_' + control + '_mgmt.txt','Fish/FishSpectra_' + control + '_mgmt')
        
    # run simulation for multiple values of h and count crashes
    h_vector = np.arange(0.01,0.11,0.01)
    if control == 'harvest':
        T2 = 100000
    else:
        T2 = 10000
    #eps = np.random.normal(0,sigma,T) # generate noise
    eps2 = np.loadtxt('./../R/epsilon2_fish.txt') # read in to verify same output from R and Python
    x = np.zeros([len(h_vector),len(u),T2])
    z = np.zeros([len(h_vector),len(u),T2])
    for i in range(len(h_vector)):
        x[i,:,:], z[i,:,:] = fishHarvest_sim(phi, h_vector[i], u, T2, eps2, 0, control, True)
        
    if control == 'noise': # run white noise simulation (phi = 0.0 and u = 0.0)
        x2 = np.zeros([len(h_vector),1,T2])
        z2 = np.zeros([len(h_vector),1,T2])
        for i in range(len(h_vector)):
            x2[i,:,:], z2[i,:,:] = fishHarvest_sim(0.0, h_vector[i], [0.0], T2, eps2, 0, control, True)
    
        x = np.concatenate((x,x2),1)
        
    # find average time till crash
    crashTime = T2/((x <= 0.12).sum(axis=2) + (x == x0).sum(axis=2) + 1)

    # plot average time till crash and compare with R
    plotCrash(n, h_vector, crashTime, labels, colors, './../R/fishCrash_' + control + '_mgmt.txt',\
        'Fish/FishCrash_' + control + '_mgmt')
    
    return None

# run simulation for harvest management
# Control parameters for harvest management
# (1 - u_harvest*B)*h*x smooths allowable catch using a weighted sum of past harvests
# -1 < u_harvest < 0  corresponds to a manager trying to minimize short-term variance
# u_harvest = 0.0 corresponds to a manager ignoring variance in their management strategy
u = np.array([-0.8, 0.0])
n = 2 # just control and no control
labels = ['Harvest Mgmt - R','Harvest Mgmt - Python','Nominal - R','Nominal - Python']
colors = ['#a50f15','#fb6a4a','#08519c','#6baed6']
runSimulation(u, n, labels, colors, 'harvest')

# run simulation for noise management
# Control parameters for noise management
# U(B) = 1 - u_noise*B smooths the noise, where U(B) is the manager's intervention and B is the backshift operator
# 0 < u_noise < 1 corresponds to a manager trying to minimize short-term variance
# u_noise = 0.0 corresponds to a manager ignoring variance in their management strategy
u = np.array([0.6, 0.0])
n = 3 # control, no control and white noise
labels = ['Noise Mgmt - R','Noise Mgmt - Python','Nominal - R','Nominal - Python', 'White Noise - R', 'White Noise - Python']
colors = ['#a50f15','#fb6a4a','#08519c','#6baed6','k','0.5']
runSimulation(u, n, labels, colors, 'noise')