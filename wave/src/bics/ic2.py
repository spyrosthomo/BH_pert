'''
    Module to define the initial conditions of the problem  
    -  u(0, xt) = sin(x*pi)
    - u_t(1,xt) =      0
'''
def icDer0():     # IC for psi 
    import numpy as np
    import inc 
    #-------------------------
    xti = inc.xti; xtf = inc.xtf;
    l   = xtf-xti
    Nxt = inc.Nxt;
    xt  = np.linspace(xti, xtf, Nxt+1)
    f1 = inc.icp1; f2 = inc.icp2; f3 = inc.icp3;
    #---
    print(l)
    pl = 2*np.sin(2.0*np.pi*f1*xt/l) + 2*np.sin(2.0*f2*np.pi*xt/l)
    return pl
#-------------------------------------------
def icDer1(ic1):
    '''
        I.C. for the derivative = 0 
        ==> u(1, xt) = u(0, xt)
    '''
    return ic1
#-------------------------------------