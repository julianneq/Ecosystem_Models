from __future__ import division
import numpy as np
from computeSpectra import computeSpectra
from make_plots import plotTimeSeries, plotSpectra, plotCrash

# Model parameters from Walker et al., 1981
# noise parameters
sigma = 5 # standard deviation of noise
phi = 0.2 # auto-regressive noise parameter

# grass (x) parameters
r_x = 4000 # growth coefficient from topsoil water
c2 = 3.15 # consumption rate by herbivore
G0 = 500 # grazing threshold
d3 = 600 # half-saturation coefficient for consumption, G
c3 = d3 - G0 #
i0 = 0.1 # minimum infiltration at low values of G
alpha = 0.122/0.21 # ratio of water uptake efficiency by wood to that of grass
L_x = 1.0 # annual grass loss by respiration and death
beta = 0.01/0.21 # ratio of rate of evaporation of soil water (as depth/time) to theta_x 
                # (what is theta_x?)
d1 = 500 # parameter determining c1, below
c1 = d1 + i0*d1/(1-i0) # rate at which infiltration increases with x
x_in = 1 # grass input (colonization) - not in original Walker model

# wood (y) parameters
r_y = 3000 # growth coefficient from topsoil water
r_S = 2000 # growth coefficient from subsoil water (4000 in Walker et al.)
d4 = 500 # half saturation effect of G on S
K = 2 # proportional reducation of r_S at high x
p = 0.1 # inefficiency of subsoil water use
L_y = 1 # annual loss of woody vegetation by respiration and death
gamma = 1/alpha #
delta = beta/alpha #

dt = 1 # length of time step

# Initial conditions
x0 = 2000 # initial grassland biomass
y0 = 2000 # initial woodland biomass

# Control parameter determining anthropogenic inputs, U(B)
# U(B)=1-u*B where U(B) is the manager's intervention and B is the backshift operator
# u = 0.07 corresponds to a manager trying to minimize short-term variance
# u = 0.0 corresponds to a manager ignoring variance in their management strategy
# u = 0.5 corresponds to a managers trying to increase short-term variance
u = np.array([0.07, 0.0, -0.5])

def grassland_sim(h, T, eps, burnin):
    # State variables
    x = np.zeros([len(u),T])
    y = np.zeros([len(u),T])
    dxdt = np.zeros([T-1])
    dydt = np.zeros([T-1])
    Q = np.zeros([T-1]) #
    z = np.zeros([T-1]) # herbivory
    
    # simulate remaining time steps for both management strategies
    # set random perturbation of initial conditions
    x[:,0] = x0 + eps[0]
    y[:,0] = y0
    for i in range(len(u)):
        for t in range(T-1):
            I = (x[i,t] + c1*i0)/(x[i,t] + c1) # infiltration
            F_x = 1/(x[i,t] + alpha*y[i,t] + beta) #
            Q[t] = max(0, c2*(1 - G0/x[i,t])/(x[i,t] + c3)) #
            S = r_S*(1 - x[i,t]/(K*(x[i,t] + d4))) #
            F_y = 1/(y[i,t] + gamma*x[i,t] + delta) #
            dxdt[t] = r_x*I*F_x*x[i,t] - L_x*x[i,t] + x_in
            dydt[t] = r_y*I*F_y*y[i,t] + S*y[i,t]/(y[i,t] + p)
            if t < 2: # just add noise
                z[t] = h*x[i,t]*Q[t]
                x[i,t+1] = x0 + eps[t+1]
                y[i,t+1] = y0
            else:
                z[t] = h*Q[t]*x[i,t] - phi*h*Q[t-1]*x[i,t-1]
                x[i,t+1] = max(1, x[i,t] + phi*(x[i,t] - x[i,t-1]) + dxdt[t]*dt - phi*dxdt[t-1]*dt \
                    - max(0.1, z[t])*dt + eps[t+1])
                y[i,t+1] = max(1, y[i,t] + dydt[t]*dt - L_y*y[i,t]*dt + u[i]*L_y*y[i,t-1])
    
    # remove burn-in
    x = x[:,burnin::]
    y = y[:,burnin::]
    
    return x, y
    
# run simulation
burnin = 500
T = 2**11 + burnin # total number of time steps
h = 10 # cattle stocking density set by rangeland manager
#eps = np.random.normal(0,sigma,T) # generate noise
eps = np.loadtxt('./../R/epsilon_grass.txt') # read in to verify same output from R and Python
x, y = grassland_sim(h, T, eps, burnin)

# plot time series of x w/ and w/o variance control and compare with R
labels = ['Noise Mgmt - R','Noise Mgmt - Python','Nominal - R','Nominal - Python', 'White Noise - R', 'White Noise - Python']
colors = ['#a50f15','#fb6a4a','#08519c','#6baed6','k','0.5']
plotTimeSeries(x, len(u), T, burnin, labels, colors, './../R/grassland_tseries.txt','Grassland/GrassTimeSeries')

# compute spectra
k = 30 # number of Slepian windows
NW = k/2 # time half bandwith parameter
Slep_real_sq, weights = computeSpectra(x, NW, k, burnin, T)

# plot spectra and compare with R
plotSpectra(len(u), T, burnin, Slep_real_sq, weights, labels, colors, './../R/grassland_spectra.txt','Grassland/GrassSpectra')

# run simulation for multiple llambdas
h_vector = np.arange(1.0, 100.0+99.0/14.0, 99.0/14.0)
T = 1500 # total number of time steps
#eps = np.random.normal(0,sigma,T) # generate noise
eps = np.loadtxt('./../R/epsilon2_grass.txt') # read in to verify same output from R and Python
x = np.zeros([len(h_vector),len(u),T-burnin])
y = np.zeros([len(h_vector),len(u),T-burnin])
for i in range(len(h_vector)):
    x[i,:,:], y[i,:,:] = grassland_sim(h_vector[i], T, eps, burnin)

# find average grass and woody vegetation as function of stocking level
x_avg = np.mean(x,2)
y_avg = np.mean(y,2)

# plot average grass and woody vegetation and compare with R
plotCrash(len(u), h_vector, x_avg, labels, colors, './../R/avg_grass.txt','Grassland/AvgGrass')
plotCrash(len(u), h_vector, y_avg, labels, colors, './../R/avg_wood.txt','Grassland/AvgWood')