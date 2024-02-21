'''
    PoschlTeller potential 
'''
#----------------------------------
class Potential:
    import inc 
    #---------------------
    def __init__(self, V0, a, x0, nulll):
        import numpy as np
        self.potentialName = "Poschl-Teller Potential"
        self.forFileNames  = "PT"
        self.saveLoc       = "../figures/PT.pdf"
        self.xti           = Potential.inc.xti
        self.xtf           = Potential.inc.xtf
        self.Nxt           = Potential.inc.Nxt
        self.V0            = V0
        self.a             = a 
        self.x0            = x0
        self.xt            = np.linspace(self.xti, self.xtf, self.Nxt + 1); 
        self.potentialForm = f"V = {V0}/cosh^2({a}(x-{x0}))"


    def potential(self):
        import numpy as np 
        #----------------------
        return self.V0/np.cosh(self.a*(self.xt-self.x0))**2
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
        ax.plot(self.xt, V, label=f'')
        ax.grid()
        plt.legend()
        #-
        fig.savefig(self.saveLoc)