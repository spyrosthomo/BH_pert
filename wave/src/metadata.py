class Metadata:
    def __init__(self, potentialName):
        import inc
        import numpy as np
        self.icp = np.array([inc.icp1, inc.icp2, inc.icp3])
        self.bcp = np.array([inc.bcp1, inc.bcp2, inc.bcp3])
        self.vp  = np.array([inc.vp1, inc.vp2, inc.vp3, inc.vp4])
        self.ti  = inc.ti;  self.tf  = inc.tf
        self.xti = inc.xti; self.xtf = inc.xtf
        self.Dt  = inc.Dt ; self.Dxt = inc.Dxt; self.Cxt = inc.Cxt
        self.lam = inc.lam
        self.potentialName = potentialName