import numpy as np
from spectrum import mtm # http://thomas-cokelaer.info/software/spectrum/html/user/install.html

def computeSpectra(x, NW, k, burnin, T):
    n = np.shape(x)[0]
    x_std = np.zeros(np.shape(x))
    Slep_complex = np.zeros([n,k,2*(T-burnin)],dtype='complex')
    Slep_real_sq = np.zeros([n,k,2*(T-burnin)])
    weights = np.zeros([n,2*(T-burnin),k])
    eigvals = np.zeros([n,k])
    for i in range(n):
        x_std[i,:] = (x[i,:]-np.mean(x[i,:]))/np.std(x[i,:])
        NFFT = int(2*2**np.ceil(np.log2(len(x_std[i,:]))))
        Slep_complex[i,:,:], weights[i,:,:], eigvals[i,:] = mtm.pmtm(x_std[i,:], NW=NW, k=k, \
            NFFT=NFFT,show=False)
        Slep_real_sq[i,:,:] = np.abs(Slep_complex[i,:,:])**2

    return Slep_real_sq, weights