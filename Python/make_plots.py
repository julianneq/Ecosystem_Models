from matplotlib import pyplot as plt
import numpy as np

def plotTimeSeries(x, n, T, burnin, labels, colors, inputFile, outputName):
    
    # load R output of time series
    R_tSeries = np.loadtxt(inputFile,skiprows=1,usecols=range(1,n+1))
    
    ax = plt.subplot(111)
    # first plot R output
    for i in range(n):
        ax.plot(R_tSeries[:,i],c=colors[i*2],label=labels[i*2],linewidth=2)

    plt.xlim([0,T-burnin-1])
    handles, legend_labels = ax.get_legend_handles_labels()
    plt.legend(handles, legend_labels,loc='upper right',ncol=2)
    plt.savefig('./../Figures/' + outputName + '_R.png')
    
    # add Python output
    for i in range(n):
        ax.plot(x[i,:],c=colors[i*2+1],label=labels[i*2+1],linewidth=2)
    
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles, labels,loc='upper right',ncol=2)
    plt.savefig('./../Figures/' + outputName + '_Rpy.png')
    plt.clf()
    
    return None
    
def plotSpectra(n, T, burnin, Slep_real_sq, weights, labels, colors, inputFile, outputName):
    
    # load R output of spectra
    var_spec_R = np.loadtxt(inputFile,skiprows=1,usecols=range(1,n+1))
    
    ax = plt.subplot(111)
    # first plot R output
    for i in range(n):
        ax.semilogy(var_spec_R[:,i],c=colors[i*2],label=labels[i*2],linewidth=2)
        
    plt.xlim([0,T-burnin-1])
    handles, legend_labels = ax.get_legend_handles_labels()
    plt.legend(handles, legend_labels, loc='lower left',ncol=2)
    plt.savefig('./../Figures/' + outputName + '_R.png')
    
    # add Python output
    for i in range(n):
        ax.semilogy(np.mean(Slep_real_sq[i,:,:]*np.transpose(weights[i,:,:]),0)[0:T-burnin],\
            c=colors[i*2+1],label=labels[i*2+1],linewidth=2)
        
    handles, legend_labels = ax.get_legend_handles_labels()
    plt.legend(handles, legend_labels, loc='lower left',ncol=2)
    plt.savefig('./../Figures/' + outputName + '_Rpy.png')
    plt.clf()

    return None
    
def plotCrash(n, x_variable, y_variable_py, labels, colors, inputFile, outputName):
    
    # load R output of crash results
    y_variable_R = np.loadtxt(inputFile,skiprows=1,usecols=range(1,n+1))
    
    ax = plt.subplot(111)
    # first plot R output
    for i in range(n):
        l1, = ax.plot(x_variable,y_variable_R[:,i],label=labels[i*2],marker='o',linewidth=2,\
            markerfacecolor='w',markeredgecolor=colors[i*2],markeredgewidth=2)
        plt.setp(l1, c=colors[i*2])
    
    handles, legend_labels = ax.get_legend_handles_labels()
    plt.legend(handles, legend_labels, loc='lower left',ncol=2)
    plt.savefig('./../Figures/' + outputName + '_R.png')
    
    # add Python output
    for i in range(n):
        l1, = ax.plot(x_variable,y_variable_py[:,i],label=labels[i*2+1],marker='o',linewidth=2,\
            markerfacecolor='w',markeredgecolor=colors[i*2+1],markeredgewidth=2)
        plt.setp(l1, c=colors[i*2+1])
        
    handles, legend_labels = ax.get_legend_handles_labels()
    plt.legend(handles, legend_labels, loc='lower left',ncol=2)
    plt.savefig('./../Figures/' + outputName + '_Rpy.png')
    plt.clf()
    
    return None