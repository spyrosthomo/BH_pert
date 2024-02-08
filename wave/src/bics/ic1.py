'''
Module to define the initial conditions of the problem  
    -  u(0, xt) = (1/sigma)*np.exp(-0.5*(xt-mean)^2/sigma^2)
    - u_t(1,xt) =      0
'''
def icDer0():        # IC for psi 
    '''
        I.C. for the solution
    '''
    import numpy as np
    import inc            #define:  (.) mean :=icp1 
                          #         (.) sigma:=icp2
    #----------------------------
    xti   = inc.xti; xtf = inc.xtf;
    Nxt   = inc.Nxt; 
    mean  = inc.icp1; 
    sigma = inc.icp2; 
    xt    = np.linspace(xti, xtf, Nxt+1)
    #--- 
    exponent = -0.5*(xt-mean)**2/sigma**2
    return (1/sigma)*np.exp(exponent); 
#-----------------------------------------------------------
def icDer1(ic1):        # IC for 1st derivative of psi 
    '''
        I.C. for the derivative = 0 
        ==> u(1, xt) = u(0, xt)
    '''
    return ic1
#---------------------------------------
def ic1Fun():
    import numpy as np 
    import inc 
    return np.zeros(inc.Nxt+1)
