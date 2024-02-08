''' 
    Module to implement Nonreflecting BCs
         UPWIND METHOD

    IN : psiJm1-> psi[j-1, :]
         psiJm2-> psi[j-2, :]
         V  -> potential
         p  -> needed parameter 

    OUT:  boundary condition value for the next time step
''' 
class BC:
    import inc
    #-----------------------
    def __init__(self, psiJm1, psiJm2=0, V=0, p=[]):
        self.psiJm1 = psiJm1
        self.psiJm2 = psiJm2
        self.V      = V
        self.p      = p
    #----------------------
    def bcF(psiJm1, psiJm2=0, V=0, p=0):
        '''
            BC for the final point
        '''
        lam = BC.inc.lam
        #--------------
        return (1-lam)*self.psiJm1[-1] + lam*self.psiJm1[-2]
    #------------------------------------------------
    def bcI(psiJm1, psiJm2=0, V=0, p=0):
        '''
            BC for the initial point
        '''
        lam = BC.inc.lam
        #--------------
        return lam*self.psiJm1[1] - (lam-1)*self.psiJm1[0]