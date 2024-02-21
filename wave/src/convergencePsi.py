import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize  import curve_fit
import sys
from scipy.fftpack import fft, ifft 
sys.path.append('./plots')
import plot3d, plot2d    
import time 
#----------------------------------------------------

Nxt     = [100, 500, 1000, 5000, 10000, 15000 ]

m       = 0
s       = 5  # ic 
sigma   = -3 # potentails
potInfo = f'GRaxial{sigma}'
icInfo  = 'Gaussian_m{}_s{}'.format(m, s) 

R       = 48

# h @ different Nxt
fig1,(ax1,ax2) = plt.subplots(2,1, sharex=True)
ax1.set_title(f'Convergence Test h @ xt={R}')
ax1.set_xlabel('t')
ax1.set_ylabel('h(t)')
ax2.set_xlabel('t')
ax2.set_ylabel('log|h(t)|')

# |h(Nxt) - h(Nxt[-1])
fig2, (ax3,ax4) = plt.subplots(2,1,sharex=True)
ax3.set_title(f'Convergence Test h @ xt={R}')
ax3.set_xlabel('t')
ax3.set_ylabel(r'dh = h(Nxt) - h(Nxt_{prev})')
ax4.set_xlabel('t')
ax4.set_ylabel('log|dh|')


for i in range(0, np.size(Nxt)): #range(0,np.size(Nxt)):
    print(f'#----------------------------------------------')
    print(f'# Nx = {Nxt[i]}')
    R       = 48
    file    = '../outputs/{}_{}_Nxt{}.npz'.format(potInfo, icInfo, Nxt[i])
    print(file)
    t1 = time.time()
    data    = np.load(file)
    psi     = data['psi']
    t       = data['t']
    xt      = data['xt']
    print(f"# Loading time: {-t1+time.time()}")
    print(f'# R before = {R}')
    NR       = np.where(xt>=R)[0][0]
    R        = xt[NR]  
    print(f'# R after = {R}')
    h        = psi[:, NR]
    logh     = np.log10(np.abs(h))

    ax1.plot(t, h   , label=f'Nx={Nxt[i]}')
    ax2.plot(t, logh, label=f'Nx={Nxt[i]}')

    if (i>0 & np.size(h)==np.size(hprev)):
        dh = h-hprev
        ax3.plot(t,               dh , label=f'{Nxt[i]}-{Nxt[i-1]}')
        ax4.plot(t,np.log(np.abs(dh)), label=f'{Nxt[i]}-{Nxt[i-1]}')
    hprev = h 

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

fig1.savefig('../figures/convergence/h_{}_{}_R{}.png'.format(potInfo, icInfo, R))
fig2.savefig('../figures/convergence/deltah_{}_{}_R{}.png'.format(potInfo, icInfo, R))



# tepak 
hmax     = np.max(h) 
maxIndex = np.where(h == hmax)
tpeak    = t[maxIndex]
print(f'tpeak={tpeak}')

# # restric t  
# tshift   = 5
# tmax     = 90
# maxIndex = np.where(t>=tmax)[0][0]
# tmin     = tpeak  + tshift
# minIndex = np.where(t>=tmin)[0][0]
# tfit     =    t[minIndex:maxIndex]
# hfit     =    h[minIndex:maxIndex]
# loghfit  = logh[minIndex:maxIndex]

plt.show()