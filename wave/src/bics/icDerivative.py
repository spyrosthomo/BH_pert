'''
    Module to define the initial conditions of the problem 
        - u(0, xt) = (1/sigma)*np.exp(-0.5*(xt-mean)^2/sigma^2)
        -u_t(0,xt) =
'''
def icDer0():
    '''
        I.C for psi
    '''
    import numpy as np 
    import inc            #define:  (.) mean :=icp1 
                          #         (.) sigma:=icp2
                          #         (.) A    :=icp3
    #----------------------------
    xti   = inc.xti; xtf = inc.xtf;
    Nxt   = inc.Nxt; 
    mean  = inc.icp1; 
    sigma = inc.icp2; 
    A     = inc.icp3;
    xt    = np.linspace(xti, xtf, Nxt+1)
    #--- 
    exponent = -0.5*(xt-mean)**2/sigma**2
    return (A/sigma)*np.exp(exponent); 
#-----------------------------------------------------------
def icDer1(ic1):        # IC for 1st derivative of u 
    '''
        I.C. for the derivative: u_t(0,x) = u_x(0,x) := g(x)
            ==> u(1, xt) = u(0, xt) + Dt*g(x) 
    '''
    import inc
    import numpy as np 
    #---------------------------------
    Dt = inc.Dt
    xti   = inc.xti; xtf = inc.xtf;
    Nxt   = inc.Nxt; 
    mean  = inc.icp1; 
    sigma = inc.icp2; 
    A     = inc.icp3;
    xt    = np.linspace(xti, xtf, Nxt+1)
    #---
    exponent = -0.5*(xt-mean)**2/sigma**2
    g        = -(A/sigma**3)*np.exp(exponent)*(xt-sigma)
    return ic1 + Dt*g