#   -   imports 
#------------------------------------------------
# -libs
import numpy as np 
import sys
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# -modules/classes etc
import inc 
sys.path.append('./bics')
import bc1, ic1                      #TODO 1.
sys.path.append('./potentials')
import harmonicOsc                   #TODO 3.
sys.path.append('./plots')
import plot3d                        #TODO 4.
sys.path.append('./methods')          
import LeapFrog                      #TODO 5.
#------------------------------------------------
print("#------------------------ Wave Equation Solver ---------------------------")
print("#-------------------------------------------------------------------------")
#-------------------------------------
#   -   inputs 
#-------------------------------------
inc.bcp1,inc.bcp2,inc.bcp3 = [float(x) for x in input("#Enter Boundary Condition Params (1-3)      : \n").split()]
inc.icp1,inc.icp2,inc.icp3 = [float(x) for x in input("#Enter Initial  Condition Params (1-3)      : \n").split()]
inc.vp1, inc.vp2, inc.vp3, inc.vp4  = [float(x) for x in input("#Enter     Potential      Params (1-4)      : \n").split()]
inc.tf , inc.xti, inc.xtf  = [float(x) for x in input("#Enter t-limit \& x*-limits: tf,  xti,  xtf : \n").split()]
inc.Nt , inc.Nxt, inc.Npl  = [int(x)   for x in input("#Enter       Nt, Nxt \& Npl                 : \n").split()]
#   all the parameters in the inc.py file are changed and ready to be used from all the other files :) 
# -easy access 
ti = inc.ti; tf  = inc.tf ; xti = inc.xti; xtf = inc.xtf; 
Nt = inc.Nt; Nxt = inc.Nxt; Npl = inc.Npl; 
#-------------------------------------

# ti  =  0 ;
# tf  =  50;
# xti = -20;
# xtf =  20;

# Nt  = 1000 
# Nxt = 100
# Nxt = int(Nt/4)
# Npl = 100

mean  = inc.icp1
sigma = inc.icp2

#   -   initializations
#-------------------------------------
# - scalars
Dt  = (tf  - ti  )/Nt ; inc.Dt  = Dt 
Dxt = (xtf - xti )/Nxt; inc.Nxt = Nxt
c   = 1 ; 
lam = c * Dt / Dxt    ; inc.lam = lam 
# - vectors
t   = np.linspace(ti , tf , Nt+1 )
xt  = np.linspace(xti, xtf, Nxt+1)
psi = np.zeros((Nt+1, Nxt+1))           # psi(t, xt)    
d2x = np.zeros(Nxt+1)
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
print("#      lambda                     = (        %.4f          )" %( inc.lam                      ))
print("#======================================================================")
#-------------------------------------
#   -   Check for convergence
#-------------------------------------
if lam > 1: 
    print("# lambda > 1: !! UNSTABLE, try other Nt, Nxt")
    sys.exit()
#-------------------------------------
#   - Apply BICs 
#-------------------------------------
#--- I.C. : u(0, xt)   = g(xt)       (1) 
#         u_t(0, xt) = h(xt)       (2)
#    We get: (.) u[0,:] from (1)
#            (.) u[1,:] from (2) as u[1,:] = u[0,:] + Dxt * icDer1    #TODO 2.
psi[0, :] = ic1.icDer0();
psi[1, :] = ic1.icDer1(psi[0, :]);                   

figIC = plt.figure()
axIC  = figIC.add_subplot(1,1,1)
axIC.set_title("IC")
axIC.plot(xt, psi[0,:])
#--- B.C.
psi[:, 0]   = 0 
psi[:, Nxt] = 0
#--- potential 
k  = 0.1
xm = (xtf-np.abs(xti))/2
V  = 0.5*k*(xt-xm)**2
V  = 0*xt
V  = harmonicOsc.potential(inc.vp1)

figPot = plt.figure()
axPot  = figPot.add_subplot(1,1,1)
axPot.set_title("0.5kx^2")
axPot.plot(xt, V)
#_-----------------------------------------
prev1 = psi[1,:]
prev2 = psi[0,:]
for j in range(2, Nt+1):
    for i in range(1,Nxt):
        d2x[i]    =  prev1[i-1] + prev1[i+1] - 2*prev1[i]
        #psi[j-1, i-1] + psi[j-1, j+1] - 2*psi[j-1, i]
    Vpart     = Dt**2*np.multiply(psi[j-1, :], V)
    psi[j, :] = 2*prev1 - prev2 + lam**2*(d2x) - Vpart
    prev2     = prev1       #psi[j-1, i]
    prev1     = psi[j, :]

#-----plots 
plot3d.plotMesh(t, xt, psi, "3d", "t", "xt", "psi")
fig2 = plt.figure()
ax2  = fig2.add_subplot(1, 1, 1)
for i in range(Nt+1):
    if i%Npl==0:
        ax2.plot(xt, psi[i,:], label=f'i={i}')
plt.legend()
#----------
# animation
fig, ax = plt.subplots()
image   = ax.imshow(psi, extent=[xt.min(), xt.max(), t.min(), t.max()], aspect='auto'
,origin='lower', interpolation='bilinear')

cbar = plt.colorbar(image)
cbar.set_label('psi(t, x)')

# Update function for animation (optional)
def update(frame):
    image.set_array(psi(t + frame * 0.1, xt))
    return image,

# Create animation (optional)
#ani = FuncAnimation(fig, update, frames=range(100), interval=50)

plt.xlabel('x')
plt.ylabel('t')
plt.title('Evolution of psi(t, x)')
plt.show()


