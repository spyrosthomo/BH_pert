'''
    class that contains everything needed for the leapfrog method in matrix formulation
'''

class Method:
    import inc 
    methodName = ""
    #-------------------------
    def __init__(self): 
        import numpy as np 
        self.methodName = "LeapFrogMat"
        #---
        lam = Method.inc.lam; 
        Nxt = Method.inc.Nxt; 
        diag     = 2.0*(1-lam**2)*np.ones(Nxt-1)
        diagNon  = lam**2*np.ones(Nxt-2)
        # diag     = 2.0*(1-lam**2)*np.ones(Nxt)
        # diagNon  = lam**2*np.ones(Nxt-1)
        diagUp   = np.insert(diagNon, 0, 0)
        diagDown = np.insert(diagNon, np.size(diagNon), 0)
        #--
        self.A   = np.eye(Nxt-1, Nxt-1, k=0)*diag + np.eye(Nxt-1, Nxt-1, k=-1)*diagDown + np.eye(Nxt-1, Nxt-1, k=1)*diagUp;
        # self.A   = np.eye(Nxt, Nxt, k=0)*diag + np.eye(Nxt, Nxt, k=-1)*diagDown + np.eye(Nxt, Nxt, k=1)*diagUp;
    #-------------------------
    def step(self, solJ, solJm1, V):
        import numpy as np 
        Dt = Method.inc.Dt
        return np.matmul(self.A, solJ[1:-1]) - solJm1[1:-1] - np.multiply(V[1:-1], solJ[1:-1])*Dt**2
        
