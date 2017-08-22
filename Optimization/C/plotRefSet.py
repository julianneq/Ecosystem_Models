import numpy as np
from matplotlib import pyplot as plt

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

# number of time steps (T) and ensemble members (S)
T = 100
S = 100

# set optimization parameters
nobjs = 2
reliability_threshold = 0.85

def lake_sim(u, llambda):
    # Initialize arrays to store objective function values and constraints
    objs = [0.0]*nobjs
    
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
                x[t+1] = max(0.1, x[t] + (u + phi)*(x[t] - x[t-1]) - u*phi*(x[t-1] - x[t-2]) \
                    + dxdt[t]*dt - (u + phi)*dxdt[t-1]*dt + u*phi*dxdt[t-2]*dt \
                    + (1 - (u + phi) + u*phi)*llambda + eps[i,t+1])
                y[t+1] = max(1.0, y[t] + dydt[t]*dt)
                if u == 0:
                    a[t-2] = llambda
                else:
                    no_u = max(0.1, x[t] + phi*(x[t] - x[t-1]) + dxdt[t]*dt \
                        - phi*dxdt[t-1]*dt + (1 - phi)*llambda + eps[i,t+1])
                    a[t-2] = llambda + (x[t+1] - no_u)
                    
                discounted_benefit[i] += alpha*a[t-2]*delta**(t-2)
                
                if x[t+1] <= 3.0: # critical P threshold
                    yrs_pcrit_met[i] += 1
            
    # compute objectives and constraints
    objs[0] = -1*np.sum(discounted_benefit)/S
    objs[1] = -1*np.sum(yrs_pcrit_met)/(T*S)
    
    return objs

# Control parameter determining anthropogenic inputs, U(B)
# U(B)=1-u*B where U(B) is the manager's intervention and B is the backshift operator
# u = 0.6 corresponds to a manager trying to minimize short-term variance
# u = 0.0 corresponds to a manager ignoring variance in their management strategy
constant_emission_objs1 = lake_sim(0.0, 0.1)
noise_mgmt_emission_objs1 = lake_sim(0.6, 0.1)

constant_emission_objs2 = lake_sim(0.0, 0.5)
noise_mgmt_emission_objs2 = lake_sim(0.6, 0.5)

constant_emission_objs3 = lake_sim(0.0, 1.5)
noise_mgmt_emission_objs3 = lake_sim(0.6, 1.5)

constant_emission_objs4 = lake_sim(0.0, 2.5)
noise_mgmt_emission_objs4 = lake_sim(0.6, 2.5)

objs = np.loadtxt('Lake_DPS.reference')


plt.scatter(-objs[:,0],-objs[:,1],c='k',s=50)
plt.scatter(-constant_emission_objs1[0],-constant_emission_objs1[1],facecolor='#08519c',edgecolor='#08519c',s=50)
plt.scatter(-constant_emission_objs2[0],-constant_emission_objs2[1],facecolor='#4292c6',edgecolor='#4292c6',s=50)
plt.scatter(-constant_emission_objs3[0],-constant_emission_objs3[1],facecolor='#41ab5d',edgecolor='#41ab5d',s=50)
plt.scatter(-constant_emission_objs4[0],-constant_emission_objs4[1],facecolor='#006d2c',edgecolor='#006d2c',s=50)

plt.scatter(-noise_mgmt_emission_objs1[0],-noise_mgmt_emission_objs1[1],facecolor='#08519c',edgecolor='#08519c',s=50,marker='s')
plt.scatter(-noise_mgmt_emission_objs2[0],-noise_mgmt_emission_objs2[1],facecolor='#4292c6',edgecolor='#4292c6',s=50,marker='s')
plt.scatter(-noise_mgmt_emission_objs3[0],-noise_mgmt_emission_objs3[1],facecolor='#41ab5d',edgecolor='#41ab5d',s=50,marker='s')
plt.scatter(-noise_mgmt_emission_objs4[0],-noise_mgmt_emission_objs4[1],facecolor='#006d2c',edgecolor='#006d2c',s=50,marker='s')

plt.xlabel('Expected Economic Benefits',fontsize=16)
plt.ylabel('Reliability',fontsize=16)
plt.tick_params(axis='both',labelsize=14)
plt.savefig('ParetoSets.png')
plt.clf()
