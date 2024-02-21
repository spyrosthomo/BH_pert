'''
    GR potential for polar perturbations of a S BH
        - l 
        - n 
        - rs = 2M => M = rs/2

'''
class Potential:
    import inc 
    #----------------------------
    def __init__(self, l, s, rs=1, nulll=0):
        import numpy as np 
        #-------
        self.potentialName = "GR Potential Polar"
        self.potentialForm = f"V = (1-rs/r)()"
        self.forFileNames  = f"GRpolar{int(n)}"
        self.saveLoc       = "../figures/GRpolar.pdf"
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
        V  = (1-rs/r)*( 2*n**2*(n+1)*r**3 + 6*n**2*(rs/2)*r**2 + 18*n*(rs/2)*r**2 + 18*n*(rs/2)**2*r \
                         + 18*(rs/2)**3 )/(n*r+ 3*rs/2 )**2/r**3 
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

