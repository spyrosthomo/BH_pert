''' 
    Module to implement Nonreflecting BCs
         UPWIND METHOD

    IN : psiJm1-> psi[j-1, :]
         psiJm2-> psi[j-2, :]
         V  -> potential
         p  -> needed parameter 

    OUT:  boundary condition value for the next time step
''' 
def bcF(psiJm1, psiJm2=0, V=0, p=[]):
    '''
        BC for the final point
   '''
    import inc 
    lam = inc.lam
    #--------------
    return (1-lam)*psiJm1[-1] + lam*psiJm1[-2]
#------------------------------------------------
def bcI(psiJm1, psiJm2=0, V=0, p=[]):
    '''
        BC for the initial point
    '''
    import inc
    lam = inc.lam
    #--------------
    return lam*psiJm1[1] - (lam-1)*psiJm1[0]