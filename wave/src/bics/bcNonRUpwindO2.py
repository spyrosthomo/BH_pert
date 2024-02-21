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
    return psiJm1[-1]-lam/2*(3*psiJm1[-1]-4*psiJm1[-2]+psiJm1[-3])
#------------------------------------------------
def bcI(psiJm1, psiJm2=0, V=0, p=[]):
    '''
        BC for the initial point
    '''
    import inc
    lam = inc.lam
    #--------------
    return psiJm1[0] + lam/2*(-psiJm1[2]+4*psiJm1[1]-3*psiJm1[0])