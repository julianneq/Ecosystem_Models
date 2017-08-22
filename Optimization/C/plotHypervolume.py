import numpy as np
from matplotlib import pyplot as plt

steps = 1000
All_HV = np.zeros([50,steps])
NFE = np.arange(200,200*steps+1,200)

for j in range(50):
    HV = np.loadtxt('./metrics/LakeDPS_S' + str(j+1) + '.metrics', skiprows=1, usecols=[0])
    if len(HV) == steps:
        All_HV[j,:] = HV
    else:
        All_HV[j,:] = np.concatenate((np.zeros(steps-len(HV)),HV),0)
    
for j in range(50):
    plt.plot(NFE, All_HV[j,:])
    
plt.xlabel('NFE',fontsize=16)
plt.ylabel('Hypervolume',fontsize=16)
plt.tick_params(axis='both',labelsize=14)
plt.savefig('HV_v_NFE.png')
plt.clf()