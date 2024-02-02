'''
    Harmonic Oscillator
'''
#----------------------------------
def potential(k):
    import numpy as np 
    import inc 
    #----------------------
    Nxt = inc.Nxt;
    xti = inc.xti; xtf = inc.xtf;
    xt  = np.linspace(xti, xtf, Nxt + 1); 
    xm  = (xtf - np.abs(xti))/2
    return 0.5*k**2*(xt[1:-1]-xm)**2
#----------------------------------    
def potentialPlot():
    import numpy as np 
    import matplotlib.pyplot as plt
    import inc 
    #----------------------
    Nxt = inc.Nxt;
    xti = inc.xti; xtf = inc.xtf;
    xt  = np.linspace(xti, xtf, Nxt + 1); 
    k   = inc.vp1; 
    V   = potential(k)
    #----------------------
    fig = plt.figure();
    ax  = fig.add_subplot(1, 1, 1)
    ax.set_xlabel("xt")
    ax.set_ylabel("V")
    ax.set_title("Harmonic Oscillator Potential, V=0.5*k*x^2")
    ax.plot(xt[1:-1], V, label=f'k={k}')
    plt.legend()
    #-
    fig.savefig("../figures/potential_plots/SHO.pdf")