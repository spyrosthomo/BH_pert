#   -   imports 
#-------------------------------------
# -libs
import numpy as np 
import importlib 
import sys
import matplotlib.pyplot as plt 
from scipy.fftpack import fft, ifft
import time
# -modules etc
import inc
import metadata
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
bcF, icF, pot, meth         = [ x for x in input("#Enter filenames for: I.C., B.C., Potential, Method:\n" ).split()]
inc.bcp1,inc.bcp2,inc.bcp3  = [float(x) for x in input("#Enter Boundary Condition Params (1-3)       : \n").split()]
inc.icp1,inc.icp2,inc.icp3  = [float(x) for x in input("#Enter Initial  Condition Params (1-3)       : \n").split()]
inc.vp1, inc.vp2, inc.vp3, inc.vp4  = [float(x) for x in input("#Enter     Potential       Params (1-4)      : \n").split()]
inc.tf , inc.xti, inc.xtf   = [float(x) for x in input("#Enter t-limit \& x*-limits: tf,  xti,  xtf  : \n").split()]
####inc.Nt , 
inc.Nxt,                    = [int(x)   for x in input("#Enter       Nt \& Nxt                       : \n").split()]
inc.Cxt,                    = [float(x) for x in input("#Enter Cxt = Dt/Dxt                          : \n").split()]
#   all the parameters in the inc.py file are changed and ready to be used from all the other files :) 
# -easy access 
ti   = inc.ti; tf  = inc.tf ; xti = inc.xti; xtf = inc.xtf; 
####Nt   = inc.Nt; 
Nxt = inc.Nxt
Cxt = inc.Cxt
icp1 = inc.icp1; icp2 = inc.icp2; icp3 = inc.icp3 
bcp1 = inc.bcp1; bcp2 = inc.bcp2; bcp3 = inc.bcp3
vp1  = inc.vp1 ; vp2  = inc.vp2 ; vp3  = inc.vp3 ; vp4 = inc.vp4
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
Dxt = (xtf - xti)/Nxt  ; inc.Dxt = Dxt
Dt  = Cxt*Dxt          ; inc.Dt = Dt
Nt  = int((tf-ti)//Dt) ; inc.Nt = Nt
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
#   - Apply ICs 
#-------------------------------------
psi[0, :]  = ic.icDer0();
psi[1, :]  = ic.icDer1(psi[0, :]);              #TODO 2.
print("# ---- Finished with I.C.")              #DEBUG
#--------------------------------------
#   - instance of the method used
#-------------------------------------
meth   = Method(); 
print("# ---- Finished with method def")                     #DEBUG
#--------------------------------------
#   - Start the iterations
#-------------------------------------
psi[1, 0  ] = bc.bcI(psi[0,:])
psi[1, Nxt] = bc.bcF(psi[0,:])

t1 = time.time()
for j in range(2, Nt+1):
    b = psi[j-1, :]
    c = psi[j-2, :]
    psi[j, 0  ]  = bc.bcI(b)
    psi[j, Nxt]  = bc.bcF(b)
    psi[j, 1:-1] = meth.step(b, c, V)
print(f"# ---- Finished with solution @ t = {time.time()-t1}")      #DEBUG

#---- save solution
potInfo = pot.forFileNames
icInfo  = ic.name()+ '_m{}_s{}'.format(int(icp1), int(icp2))
### com     = f"# Psi for  *IC : {icInfo} : icp1={icp1}, icp2={icp2}, icp3={icp3} \n \
###             #    \t       Pot: {potInfo}: vp1 ={vp1} , vp2 ={vp1} , vp3 ={vp3}, vp4 ={vp4} \n \
###             #    \t       xti = {xti} \t xtf = {xtf} \n \
###             #    \t       tf  = {tf}"
meta    = metadata.Metadata(potInfo)
t1      = time.time()
np.savez('../outputs/{}_{}_Nxt{}.npz'.format(potInfo, icInfo, Nxt), psi=psi, t=t, xt=xt, metadata=meta)
print(f'# ---- Finished with save @ t = {time.time()-t1}')
#--------------------------------------
#   - Plots
#-------------------------------------
t1 = time.time()
# 1. Potential
pot.potentialPlot(V)

# 2. <- killed

# 3. plot IC                                 #TODO 4.
plot2d.plot2d(xt, psi[0, :], "Initial Condition", 'xt', "\psi(t=0,x)", "", "../figures/IC.pdf")

# 4. 2d plots                                #TODO 4.
fig2 = plt.figure()
ax2  = fig2.add_subplot(1, 1, 1)
for i in range(Nt+1):
    if i%((Nt+1)//5) == 0 :
        ax2.plot(xt, psi[i, :], label=f'i={i}')
plt.legend()
ax2.set_xlabel('xt')
ax2.set_ylabel('psi(n*t, xt)')
plot3d.saveFig(fig2, "../figures/2d.pdf")

#------------
# 4.5 : x - xt 
if pot.potentialName == "GR Potential" :
    pot.plotrsr()

# 5. colors
fig, ax = plt.subplots()
image   = ax.imshow(psi, extent=[xt.min(), xt.max(), t.min(), t.max()], aspect='auto'
        ,origin='lower', interpolation='bilinear')
cbar = plt.colorbar(image)
cbar.set_label('psi(t, x)')
ax.set_xlabel('xt')
ax.set_ylabel('t')
plot2d.saveFig(fig, '../figures/2d_color.pdf')

#7. xm = (xtf-xti)/2
#   plot   psi(:, xm) 
###fig = plt.figure()
###ax1 = fig.add_subplot(2,1,1)
###ax2 = fig.add_subplot(2,1,2)
#-------
##Nxtm    = int(5*Nxt//(8)) 
##xtm     = xt[Nxtm]#(xtf - np.abs(xti))/2
###xtObserverIndex = np.where(xt>=xtObserver)[0][0]
###psiM    = psi[:, xtObserverIndex]#psi[:, Nxtm]
###saveLoc = '../figures/tEvolFFT.pdf'
#-------
###ax1.plot(t, (psiM), 'r')
###ax1.set_title(f'Time evolution, @x={xtObserver}')
###ax1.set_xlabel('t')
###ax1.set_ylabel("\psi(t, x=(xf-xi)/2)")
###ax2.set_xlabel("omega(Hz)")
###ax2.set_ylabel("FFT Ampl |X(omega)|")

#8. FFT of 7.
###N     = np.size(psiM)
###tstep = tf/Nt
###Fs    = 1/tstep
###tt    = np.linspace(ti, tf         , N)
###f     = np.linspace(0 , (N-1)*tstep, N)
###Y     = fft(psiM)
###Y_mag = np.abs(Y) / N

###Ns            = N//2+1
###f_plot        = np.linspace(0, Fs/2, Ns)
###df            = 0.5*Fs/(Ns)
###y_mag_plot    = 2*Y_mag[0:Ns]
###y_mag_plot[0] = y_mag_plot[0]/2
# - print freqs
###lowerLim = 0.2
###indices  = np.where(y_mag_plot>lowerLim)[0]
###amps     = y_mag_plot[indices]
###freqs    = indices*df

#----
#ax2.stem(f_plot, y_mag_plot, 'b', markerfmt='', basefmt='-b')
#ax2.legend()
#------------
# save
#   -fft
###output1 = np.hstack((f_plot.reshape(Ns,1), y_mag_plot.reshape(Ns,1))) #np.vstack((amps[:,np.newaxis], freqs[:,np.newaxis]))
###np.savetxt(f'../outputs/fft_{pot.forFileNames}R{int(xtObserver)}.txt', output1)
###plot2d.saveFig(fig, '../figures/fft.pdf')
#   -time evol
###output2 = np.hstack((t.reshape(Nt+1,1), psiM.reshape(Nt+1,1) ))
###np.savetxt(f'../outputs/psiM{pot.forFileNames}R{int(xtObserver)}.txt', output2)

#--- indicative plot of h(t) for two observers @ 3L/8 and 5L/8
fig = plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
#-------
Nxtm1       = int(3*Nxt//(8)) 
xtObs1      = round(xt[Nxtm1],0)
xtObsIndex1 = np.where(xt>=xtObs1)[0][0]
psiObs1     = psi[:, xtObsIndex1]
Nxtm2       = int(5*Nxt//(8)) 
xtObs2      = round(xt[Nxtm2])
xtObsIndex2 = np.where(xt>=xtObs2)[0][0]
psiObs2     = psi[:, xtObsIndex2]
saveLoc     = '../figures/h_2obs.pdf'
#-------
ax1.plot(t, psiObs1, '-r')
ax2.plot(t, psiObs2, '-b')
ax1.set_title(f'Time evolution, @x={xtObs1}')
ax1.set_xlabel('t')
ax1.set_ylabel(f"\psi(t, xt={xtObs1})")
ax2.set_title(f'Time evolution, @x={xtObs2}')
ax2.set_xlabel('t')
ax2.set_ylabel(f"\psi(t, xt={xtObs2})")

#------------ TEST BOUNDARY NONREFLECTING 
plot2d.plot2d(t, psi[:, 0], f"BC @ x={xti}", "t", f"psi(t,{xti})", '', '../figures/BCleft.pdf')
plot2d.plot2d(t, psi[:,-1], f"BC @ x={xtf}", "t", f"psi(t,{xtf})", '', '../figures/BCright.pdf')

#---------------------
#plt.show()
print(f'# ---- Finished with plots @ t = {time.time()-t1}')