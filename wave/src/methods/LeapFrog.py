'''
    class that implements the leapfrog mehod
''' 
class Method:
    import inc
    methodName = ""
    #-----------------
    def __init__(self):
        self.methodName = "LeapFrog"
    #-----------------
    def step(self, prev1, prev2, V):
        import numpy as np 
        #---
        Nxt = Method.inc.Nxt;
        Dt  = Method.inc.Dt;
        lam = Method.inc.lam
        d2x = np.zeros(Nxt+1)
        # compute x-derivative
        for i in range(1, Nxt):
            d2x[i] = prev1[i-1] + prev1[i+1] - 2*prev1[i]
        # the rest
        Vpart      = Dt**2*np.multiply(prev1, V)
        return 2*prev1 - prev2 + lam**2*(d2x) - Vpart