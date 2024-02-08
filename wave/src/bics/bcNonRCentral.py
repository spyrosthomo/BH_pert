'''
    Module to implement the BCs of Nonreflecting boundaries
'''
def bcF(ic0=0, ic1=0, V=0, p=[]):
    '''
        BC for the final point xtf 
    '''
    import numpy as np 
    import inc 
    #------------------------
    Nt  = inc.Nt 
    Dt  = inc.Dt
    lam = inc.lam
    bc  = np.zeros(Nt+1)
    #-----
    # COMPLICATED STUFF FOR NO REASON ( -|cosfonor|- ), START:
    #
    # --to compute bc[1]
    # a       = (2-lam)/(lam+1)
    # b       = 2*(1-lam)
    # c       = 2*lam**2/(1+lam)
    # d       = (1-lam)/(1+lam)
    # e       = Dt**2/(lam+1)
    # bc0Left = ic0[Nt-1]           # . . . o   j=1
    #                               # . . . o   j=0               => I.C.
    #                               #   ^             : bc0Left 
    #^^^^^^^^---- END of cosfonor ----  

    # --first 2 by hand 
    bc[0] = ic0[-1]
    bc[1] = ic1[-1]           #b/c*bc[0] + c/a*bc0Left + d/a*Dt*
    
    # --all the others iteratively
    # ---some constants needed
    a = 2*(1-lam)
    b = 2*lam**2/(1+lam)
    c = (lam-1) /(lam+1)
    d =  Dt**2  /(lam+1)
    for j in range(2, Nt+1):
        bcMLeft = 

#-----------------------------------
def bcI(ic0=0, ic1=0, V=0, p=[]):
    ''' 
        BC for the initial point xti
    '''
    import numpy as np 
    import inc 
    #-------------------------
    Nt = inc.Nt 
