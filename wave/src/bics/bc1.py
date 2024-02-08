'''
    Module to define the boundary condition of the problem:
        - u(t, 0)   = 0
        - u(t, Nxt) = 0 
'''
def bcI(ic0=0, ic1=0, V=0):
    '''
        B.C. for the initial point
    '''
    import numpy as np 
    import inc
    #----------------------
    Nt = inc.Nt
    return np.zeros((Nt+1))
#---------------------------------------
def bcF(ic0=0, ic1=0, V=0):
    '''
        B.C. for the final point
    '''
    import numpy as np 
    import inc
    #----------------------
    Nt = inc.Nt
    return np.zeros((Nt+1))