import numpy as np
import matplotlib.pyplot as plt

refSet = np.loadtxt('FishGame.reference')
plt.scatter(refSet[:,0], refSet[:,1])
plt.xlabel('Negative NPV of Harvest')
plt.xlim((1.1*np.min(refSet[:,0],0)))
plt.ylabel('St Dev of Harvest')
plt.ylim((0,1.1*np.max(refSet[:,1])))
plt.savefig('FishGame_RefSet.png')