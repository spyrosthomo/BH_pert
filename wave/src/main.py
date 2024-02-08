#   -   imports 
#-------------------------------------
# -libs
import numpy as np 
import importlib 
import sys
import matplotlib.pyplot as plt 
from scipy.fftpack import fft, ifft
# from numpy.fft import fft, ifft, fftfreq
# -modules etc
import inc
sys.path.append('./plots')
import plot3d, plot2d    
sys.path.append('./bics') 
sys.path.append('./potentials')
sys.path.append('./methods')        
# ----> finish imports after the reading of the arguments
print("#------------------------ Wave Equation Solver ---------------------------")
print("#-------------------------------------------------------------------------")
#-------------------------------------
#   -   inputs 
#-------------------------------------
bcF, icF, pot, meth           = [ x for x in input("#Enter filenames for: I.C., B.C., Potential, Method:\n" ).split()]
inc.bcp1,inc.bcp2,inc.bcp3  = [float(x) for x in input("#Enter Boundary Condition Params (1-3)       : \n").split()]
inc.icp1,inc.icp2,inc.icp3  = [float(x) for x in input("#Enter Initial  Condition Params (1-3)       : \n").split()]
inc.vp1, inc.vp2, inc.vp3, inc.vp4  = [float(x) for x in input("#Enter     Potential       Params (1-4)      : \n").split()]
inc.tf , inc.xti, inc.xtf   = [float(x) for x in input("#Enter t-limit \& x*-limits: tf,  xti,  xtf  : \n").split()]
inc.Nt , inc.Nxt            = [int(x)   for x in input("#Enter       Nt \& Nxt                       : \n").split()]
#   all the parameters in the inc.py file are changed and ready to be used from all the other files :) 
# -easy access 
ti      = inc.ti; tf  = inc.tf ; xti = inc.xti; xtf = inc.xtf; 
Nt      = inc.Nt; Nxt = inc.Nxt; Npl = inc.Npl; 
# ---- Finish imports 
ic      = importlib.import_module(icF)
bc      = importlib.import_module(bcF)
#    -pot
fpath     = './potentials/{}.py'.format(pot)
spec      = importlib.util.spec_from_file_location(pot, fpath)
module    = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
Potential = getattr(module, 'Potential')
#    -meth
fpath     = './methods/{}.py'.format(meth)
spec      = importlib.util.spec_from_file_location(meth, fpath)
module    = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
Method    = getattr(module, 'Method') 
#-------------------------------------
#   -   initializations
#-------------------------------------
# - scalars
Dt     = (tf  - ti )/Nt   ; inc.Dt  = Dt 
Dxt    = (xtf - xti)/Nxt ; inc.Dxt = Dxt
#Dt     = 0.25*Dxt
c      = 1.0 ;                 # :) 
lam    =  c * Dt / Dxt   ; inc.lam = lam
# - vectors 
t      = np.linspace(ti , tf , Nt  +1)
xt     = np.linspace(xti, xtf, Nxt +1)
psi    = np.zeros((Nt+1,Nxt+1))                 # psi(t, xt)
#-------------------------------------
#   -   print all the useless 
#-------------------------------------
print("#======================================================================")
print("# Files BICS      :  ic,  bc -> %s, %s " %(icF, bcF))
print("# Potential       :    %s"               %(pot))
print("# Method          :    %s"               %(meth))
print("# Time  interval  : (ti , tf )    = (%.3f, %.3f) " %(ti , tf ))
print("# Space interval  : (xti, xtf)    = (%.3f, %.3f) " %(xti, xtf))
print("# Steps           : Nt = %d,  Nxt = %d      " %(Nt, Nxt ))
print("# potential params:     vp1-4     = (%.4f, %.4f, %.4f, %.4f)" %(inc.vp1, inc.vp2, inc.vp3, inc.vp4))
print("#   I.C.    params: icp1-3        = (   %.4f, %.4f, %.4f   )" %( inc.icp1, inc.icp2, inc.icp3 ))
print("#   B.C.    params: bcp1-3        = (   %.4f, %.4f, %.4f   )" %( inc.bcp1, inc.bcp2, inc.bcp3 ))
print("#      lambda                     = (        %.4f          )" %( inc.lam                      ))
print("#=======================================================================")
#-------------------------------------
#   -   Check for convergence
#-------------------------------------
if lam > 1: 
    print("# lambda > 1: !! UNSTABLE, try other Nt, Nxt")
    sys.exit()
#--------------------------------------
#   - Compute the potential vector
#--------------------------------------
pot = Potential(inc.vp1, inc.vp2, inc.vp3, inc.vp4)
V   = pot.potential()
print("# ---- Finished with V definition")   #DEBUG
#-------------------------------------
#   - Apply BICs 
#-------------------------------------
# - I.C.: u(0, xt)   = g(xt)       (1) 
#         u_t(0, xt) = h(xt)       (2)
#    We get: (.) u[0,:] from (1)
#            (.) u[1,:] from (2) as u[1,:] = u[0,:] + Dxt * icDer1
psi[0, :]  = ic.icDer0();
psi[1, :]  = ic.icDer1(psi[0, :]);                       #TODO 2.
print("# ---- Finished with I.C.")              #DEBUG
# - B.C. 
#psi[:,   0] = bc.bcI();             
#psi[:, Nxt] = bc.bcF();
print("# ---- Finished with B.C.")              #DEBUG
#--------------------------------------
#   - instance of the method used
#-------------------------------------
meth   = Method(); 
print("# ---- Finished with method def")                     #DEBUG
#--------------------------------------
#   - Start the iterations
#-------------------------------------
psi[1,  0] = 2*(1-lam**2)/(2-lam)*psi[0,0]  + 2*lam**2/(2-lam)*psi[0, 1] - (lam-1)/(lam+1)*2*Dt*0            - Dt**2/(2-lam)*V[0 ]*psi[0, 0]
psi[1,Nxt] = 2*(1-lam**2)/(2-lam)*psi[0,-1] + 2*lam**2/(2-lam)*psi[0,-2] + (1-lam**2)/(2-lam)/(2-lam)*2*Dt*0 - Dt**2/(2-lam)*V[-1]*psi[0,-1]
#psi[1, 0  ] = bc.bcI(psi[0,:])#lam*psi[0,1]      - (lam-1)*psi[0, 0]
#psi[1, Nxt] = bc.bcF(psi[0,:])#(1-lam)*psi[0,-1] + lam*psi[0, -2]
for j in range(2, Nt+1):
    b = psi[j-1, :]
    c = psi[j-2, :]
    # proxeirh implementation of NonR   ---------- cosfonor START |
    # right
#    psi[j, Nxt] = (2*(1-lam)*b[-1] + 2*lam**2/(1+lam)*b[-2] + (lam-1)/(lam+1)*c[-1] -
#        Dt**2/(lam+1)*V[-1]*b[-1])
    # left
#    psi[j,  0]  = (2*(1-lam)*b[0] + 2*lam**2/(1+lam)*b[1] - (lam-1)/(1+lam)*c[0] -
#        Dt**2/(1+lam)*V[0]*b[0])
    #  -------------------------------------------- cosfonor END |^^^^
    psi[j, 0  ] = bc.bcI(b)#lam*b[1]      - (lam-1)*b[0]
    psi[j, Nxt] = bc.bcF(b)#(1-lam)*b[-1] +     lam*b[-2]
    if (meth.methodName   == "LeapFrog"):
        psi[j, 1:-1]    = meth.step(b, c, V)[1:-1]
    elif (meth.methodName == "LeapFrogMat"):
        psi[j, 1:-1] = meth.step(b, c, V)

print("# ---- Finished with solution")      #DEBUG
#--------------------------------------
#   - Plots
#-------------------------------------
# 1. Potential
pot.potentialPlot(V)
# 2. 3d 
saveLocMesh  = "../figures/grid3d.pdf"
title        = "Wave equation"
xlabel       = "t"
ylabel       = "xt"
zlabel       = "psi"
plot3d.plotMesh(t, xt, psi, title, xlabel, ylabel, zlabel, 
        saveLocMesh, show=0)
# 3. plot IC                           #TODO 4.
plot2d.plot2d(xt, psi[0, :], "Initial Condition", 'xt', "\psi(t=0,x)", "", "../figures/IC.pdf")

# 4. 2d plots                                #TODO 4.
fig2 = plt.figure()
ax2  = fig2.add_subplot(1, 1, 1)
for i in range(Nt+1):
    if i%((Nt+1)//10) == 0 :
        ax2.plot(xt, psi[i, :], label=f'i={i}')
plt.legend()
ax2.set_xlabel('xt')
ax2.set_ylabel('psi(n*t, xt)')
plot3d.saveFig(fig2, "../figures/2d.pdf")

# 5. colors
fig, ax = plt.subplots()
image   = ax.imshow(psi, extent=[xt.min(), xt.max(), t.min(), t.max()], aspect='auto'
        ,origin='lower', interpolation='bilinear')
cbar = plt.colorbar(image)
cbar.set_label('psi(t, x)')
ax.set_xlabel('xt')
ax.set_ylabel('t')
plot2d.saveFig(fig, '../figures/2d_color.pdf')
#==============================================

# 6. energy: t-E(t)
                            #TODO 6.

#7. xm = (xtf-xti)/2
#   plot   psi(:, xm) 
fig = plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
#-------
xtm     = (xtf - np.abs(xti))/2
Nxtm    = Nxt//2 + 1 
psiM    = psi[:, Nxtm]#psi[:, Nxtm]
saveLoc = '../figures/tEvolFftMid.pdf'
np.savetxt('../outputs/psiM.txt', psiM)
#-------
ax1.plot(t, psiM, 'r')
ax1.set_title(f'Time evolution in the middle of the interval, @x={xtm}')
ax2.set_xlabel('t')
ax2.set_ylabel("\psi(t, x=(xf-xi)/2)")

#8. FFT of 7.
N     = np.size(psiM)
tstep = tf/Nt
Fs    = 1/tstep

tt    = np.linspace(ti, tf         , N)
f     = np.linspace(0 , (N-1)*tstep, N)

Y     = fft(psiM)
Y_mag = np.abs(Y) / N

f_plot        = np.linspace(0, Fs/2, N//2+1)
df            = 0.5*Fs/(N//2+1)
y_mag_plot    = 2*Y_mag[0:N//2+1]
y_mag_plot[0] = y_mag_plot[0]/2
# - print freqs
lowerLim = 0.2
indices  = np.where(y_mag_plot>lowerLim)[0]
amps     = y_mag_plot[indices]
freqs    = indices*df
np.savetxt('../outputs/fft_output.txt', np.hstack((amps[:,np.newaxis], freqs[:,np.newaxis])))

#f_mes    = f'amp={(y_mag_plot[indices])}, freq={indices*df}'
#print(f_mes)
#----

ax1.plot(t, psiM, 'r')
ax1.set_title(f'Time evolution in the middle of the interval, @x={xtm}')
ax2.set_xlabel('t')
ax2.set_ylabel("\psi(t, x=(xf-xi)/2)")
#----

ax2.stem(f_plot, y_mag_plot, 'b', markerfmt='', basefmt='-b')
ax2.set_xlabel("Freq(Hz)")
ax2.set_ylabel("FFT Ampl |X(Freq)|")
ax2.legend()
#--------
plot2d.saveFig(fig, '../figures/fft.pdf')


#------------ TEST BOUNDARY NONREFLECTING 
plot2d.plot2d(t, psi[:, 0], f"BC @ x={xti}", "t", f"psi(t,{xti})", '', '../figures/BCleft.pdf')
plot2d.plot2d(t, psi[:,-1], f"BC @ x={xtf}", "t", f"psi(t,{xtf})", '', '../figures/BCright.pdf')

#---------------------
plt.show()