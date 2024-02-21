'''
    Module for 2d plotting
'''
def plot2d(x, y, title, xlabel, ylabel, leg='', saveLocation=0, show=0):
    '''
        simple 2d plots
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    #----------------------------
    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)
    ax.plot(x, y, label=leg)
    # -- cosmetics
    ax.set_title(title)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.grid()
    # -- save 
    if (saveLocation):
        saveFig(fig, saveLocation)
    if show:
        plt.show(block=1)

#--------------------------------
def saveFig(fig, saveLocation): 
    ''' 
        Saves the figure
    '''
    import matplotlib.pyplot as plt
    fig.savefig(saveLocation)        