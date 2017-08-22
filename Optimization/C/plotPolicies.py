import numpy as np
from matplotlib import pyplot as plt

objs = np.loadtxt('Lake_DPS.reference')
resultfile = np.loadtxt('Lake_DPS.resultfile')

bestEconSoln = np.argmin(objs[:,0])
bestRelSoln = np.argmin(objs[:,1])
medianRelSoln = np.argsort(objs[:,1])[int(len(objs[:,1])/2)]

vars_econ = resultfile[bestEconSoln,0:4]
vars_rel = resultfile[bestRelSoln,0:4]
vars_med = resultfile[medianRelSoln,0:4]

x = np.arange(0.0,8.1,0.1)
a_econ = np.zeros(len(x))
a_rel = np.zeros(len(x))
a_med = np.zeros(len(x))
for i in range(len(x)):
    a_econ[i] = min(1.5, max(max(vars_econ[0]*x[i] + vars_econ[1], vars_econ[2]*x[i] + vars_econ[3]), 0.01))
    a_rel[i] = min(1.5, max(max(vars_rel[0]*x[i] + vars_rel[1], vars_rel[2]*x[i] + vars_rel[3]), 0.01))
    a_med[i] = min(1.5, max(max(vars_med[0]*x[i] + vars_med[1], vars_med[2]*x[i] + vars_med[3]), 0.01))
    
plt.plot(x,a_econ,c='#006d2c',linewidth=2)
plt.plot(x,a_rel,c='#08519c',linewidth=2)
plt.plot(x,a_med,c='k',linewidth=2)
plt.plot([3,3],[0,3],c='#e41a1c',linewidth=2)
plt.xlabel('Lake P Concentration',fontsize=16)
plt.ylabel('P Emission',fontsize=16)
plt.tick_params(axis='both',labelsize=14)
plt.savefig('Policies.png')
plt.clf()