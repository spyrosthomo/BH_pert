'''
    Harmonic Oscillator
'''
#----------------------------------
class Potential:
    import inc 
    #---------------------
    def __init__(self, k):
        import numpy as np
        self.potentialName = "Harmonic Oscillator Potential"
        self.potentialForm = "V = 0.5*k*x^2"
        self.saveLoc       = "../figures/SHO.pdf"
        self.k             = k
        self.xti           = Potential.inc.xti
        self.xtf           = Potential.inc.xtf
        self.Nxt           = Potential.inc.Nxt
        self.xt            = np.linspace(self.xti, self.xtf, self.Nxt + 1); 


    def potential(self):
        import numpy as np 
        #----------------------
        # Nxt = inc.Nxt;
        # xti = inc.xti; xtf = inc.xtf;
        xm  = (self.xtf - np.abs(self.xti))/2
        return 0.5*self.k**2*(self.xt-xm)**2
#----------------------------------    
    def potentialPlot(self, V):
        import numpy as np 
        import matplotlib.pyplot as plt
        #----------------------
        fig = plt.figure();
        ax  = fig.add_subplot(1, 1, 1)
        ax.set_xlabel("xt")
        ax.set_ylabel("V")
        ax.set_title(self.potentialName + ", " + self.potentialForm)
        ax.plot(self.xt, V, label=f'k={self.k}')
        ax.grid()
        plt.legend()
        #-
        fig.savefig(self.saveLoc)