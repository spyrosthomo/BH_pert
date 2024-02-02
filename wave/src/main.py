#   -   imports 
#-------------------------------------
# -libs
import numpy as np 
import importlib 
import sys
# -modules etc
import inc
sys.path.append('./bics') 
import bc1, ic1                      #TODO 1.
sys.path.append('./potentials')
sys.path.append('./plots')
import plot3d                        #TODO 3.

#-------------------------------------
print("#------------------------ Wave Equation Solver ---------------------------")
print("#-------------------------------------------------------------------------")
#-------------------------------------
#   -   inputs 
#-------------------------------------
inc.bcp1,inc.bcp2,inc.bcp3 = [float(x) for x in input("#Enter Boundary Condition Params (1-3)      : \n").split()]
inc.icp1,inc.icp2,inc.icp3 = [float(x) for x in input("#Enter Initial  Condition Params (1-3)      : \n").split()]
inc.vp1, inc.vp2, inc.vp3, inc.vp4  = [float(x) for x in input("#Enter     Potential      Params (1-4)      : \n").split()]
inc.tf , inc.xti, inc.xtf  = [float(x) for x in input("#Enter t-limit \& x*-limits: tf,  xti,  xtf : \n").split()]
inc.Nt , inc.Nxt           = [int(x)   for x in input("#Enter          Nt \& Nxt                   : \n").split()]
#   all the parameters in the inc.py file are changed and ready to be used from all the other files :) 
# -easy access 
ti = inc.ti; tf = inc.tf; xti = inc.xti; xtf = inc.xtf; 
Nt = inc.Nt; Nxt = inc.Nxt;
#-------------------------------------
#   -   initializations
#-------------------------------------
# - scalars
Dt     = (tf  - ti)/Nt   ; inc.Dt  = Dt 
Dxt    = (xtf - xti)/Nxt ; inc.Dxt = Dxt
c      = 1 ;                 # :) 
lam    =  c * Dt / Dxt ;
# - vectors 
t      = np.linspace(0  , tf , Nt  +1)
xt     = np.linspace(xti, xtf, Nxt +1)
psi    = np.zeros((Nt+1,Nxt+1))                 # psi(t, xt)
#-------------------------------------
#   -   print all the useless 
#-------------------------------------
print("#======================================================================")
print("# Time  interval  : (ti , tf )    = (%.3f, %.3f) " %(ti , tf ))
print("# Space interval  : (xti, xtf)    = (%.3f, %.3f) " %(xti, xtf))
print("# Steps           : Nt = %d,  Nxt = %d      " %(Nt, Nxt ))
print("# potential params:     vp1-4     = (%.4f, %.4f, %.4f, %.4f)" %(inc.vp1, inc.vp2, inc.vp3, inc.vp4))
print("#   I.C.    params: icp1-3        = (   %.4f, %.4f, %.4f   )" %( inc.icp1, inc.icp2, inc.icp3 ))
print("#   B.C.    params: bcp1-3        = (   %.4f, %.4f, %.4f   )" %( inc.icp1, inc.icp2, inc.icp3 ))
print("#======================================================================")
#-------------------------------------
#   - Apply BICs 
#-------------------------------------
# - I.C.: u(0, xt)   = g(xt)       (1) 
#         u_t(0, xt) = h(xt)       (2)
#    We get: (.) u[0,:] from (1)
#            (.) u[1,:] from (2) as u[1,:] = u[0,:] + Dxt * icDer1
psi[0, :] = ic1.icDer0();
psi[1, :] = ic1.icDer1();                       #TODO 2.
print("# ---- Finished with I.C.")              #DEBUG
# - B.C. 
print(np.shape(psi))                            #DEBUG
psi[:,   0  ] = bc1.bcI();             
psi[:, Nxt] = bc1.bcF();
print("# ---- Finished with B.C.")              #DEBUG
#--------------------------------------
#   - Declare/Define the constant matrix 
#-------------------------------------

diag     = 2.0*(1-lam**2)*np.ones(Nxt-1)
diagNon  = lam**2*np.ones(Nxt-2)
diagUp   = np.insert(diagNon, 0, 0)
diagDown = np.insert(diagNon, np.size(diagNon), 0)

#--
A = np.eye(Nxt-1, Nxt-1, k=0)*diag + np.eye(Nxt-1, Nxt-1, k=-1)*diagDown + np.eye(Nxt-1, Nxt-1, k=1)*diagUp;

print("# ---- Finished with A def")                     #DEBUG
#--------------------------------------
#   - Compute the potential vector
#-------------------------------------
V = np.zeros((Nxt-1))
print("# ---- Finished with V definition")   #DEBUG
#--------------------------------------
#   - Start the iterations
#-------------------------------------
for j in range(1, Nt):
    b = psi[j  , 1:-1];
    c = psi[j-1, 1:-1];
    psi[j+1, 1:-1] = np.matmul(A, b) - c - V


print("# ---- Finished with solution")      #DEBUG
# print(np.shape(t))              #DEBUG  
# print(np.shape(xt))             #DEBUG
# print(np.shape(psi))            #DEBUG

#--------------------------------------
#   - Plot
#-------------------------------------
saveLocMesh  = "../figures/grid3d.pdf"
title        = "Wave equation"
xlabel       = "t"
ylabel       = "xt"
zlabel       = "psi"
plot3d.plotMesh(t, xt, psi, title, xlabel, ylabel, zlabel, 
        saveLocMesh, save=1, show=0)

# 2d-plot                                #TODO 4.
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
#-----------------------------------------------
fig2 = plt.figure()
ax2  = fig2.add_subplot(1, 1, 1)
for i in range(Nt+1):
    if i%10 == 0 :
        ax2.plot(xt, psi[i, :], label=f'i={i}')
plt.legend()
plt.show(block=True)

plot3d.saveFig(fig2, "../figures/2d.pdf")