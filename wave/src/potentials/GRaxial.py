'''
    GR potential for axial perturbations of a S BH
        - l 
        - sigma       : compactification of 1.) sigma = 0: electromagnetic perturbations
                                            2.) sigma = 1: scalar          perturbations
                                            3.) sigma =-3: gravitational   perturbations [kokkotas]
                            sigma=1-spin^2
        - rs: Horizon
'''
class Potential:
    import inc 
    #----------------------------
    def __init__(self, l, s, rs=1, nulll=0):
        import numpy as np 
        #-------
        self.potentialName = "GR Potential"
        self.potentialForm = f"V = (1-rs/r)({l}({l}+1)/r^2 + {s}rs/r^3)"
        self.forFileNames  = f"GRaxial{int(s)}"
        self.saveLoc       = "../figures/GRtestField.pdf"
        self.xti           = Potential.inc.xti
        self.xtf           = Potential.inc.xtf
        self.Nxt           = Potential.inc.Nxt
        self.l             = l 
        self.rs            = rs
        self.s             = s
        self.xt            = np.linspace(self.xti, self.xtf, self.Nxt+1)
    #-----------------------------
    def invert(self):
        ''' 
            For given r* find r 
                - X:=exp((r-rs)/rs)
                - Y:=(r-rs)/rs
        '''
        import scipy 
        import numpy as np
        X = np.exp((self.xt-self.rs)/self.rs)
        Y = np.zeros(self.Nxt+1)
        for i in range(np.size(X)): 
            if   (X[i] >= 0):          
                Y[i] = scipy.special.lambertw(X[i], 0, tol=1e-8)
            elif (X[i] <  0): 
                Y[i] = scipy.special,lambertw(X[i],-1, tol=1e-8)
        x = self.rs*Y + self.rs
        return x
    #----------------------------
    def potential(self):
        import numpy as np 
        rs = self.rs
        l  = self.l
        r  = Potential.invert(self)
        V  = (1-rs/r)*( l*(l+1) + self.s*rs/r)/r**2 
        return V 
    #----------------------------
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
    #------------------------------
    def plotrsr(self, save=0):
        import matplotlib as plt 
        import plot2d
        #----------------------------------
        xt = self.xt
        x  = self.invert()
        plot2d.plot2d(x, xt, "xtort (x)", "x", "x_t", '', saveLocation='../figures/rVSrs.pdf')

