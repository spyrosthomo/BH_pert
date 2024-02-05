'''
    Module for 3d plotting 
'''
def plotMesh(x, y, u, title, xlabel, ylabel, zlabel, saveLocation=0, show=0):
    '''
        plots 3d surfs on a grid 
    '''
    import numpy as np 
    import matplotlib.pyplot as plt 
    from mpl_toolkits.mplot3d import Axes3D
    #-----------------------------------------------
    X ,Y = np.meshgrid(x, y)
    X = np.transpose(X); Y = np.transpose(Y)
    fig  = plt.figure();
    ax   = fig.add_subplot(projection='3d')
    ax.plot_surface(X, Y, u)
    #   -- cosmetics
    ax.set_title(title)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_zlabel(zlabel);
    #   -- save
    if saveLocation:
        saveFig(fig, saveLocation)
    if show:
        plt.show(block=1)
#---------------------------------------------------
def saveFig(fig, saveLocation):
    '''
        Saves the figure
    '''
    import matplotlib.pyplot as plt
    fig.savefig(saveLocation)