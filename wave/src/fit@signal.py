import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize  import curve_fit
import sys
from scipy.fftpack import fft, ifft 
sys.path.append('./plots')
import plot3d, plot2d    
import time 
#--------------------------------------
def modelQNMorg(t, A, oR, oI, phi):
    p = A*np.exp(-oI*(t))*np.sin(oR*(t)+phi) 
    return p
#-----------------------------------------------
def modelQNM2(t, A1, oR1, oI1, phi1, A2, oR2, oI2, phi2):
    p = A1*np.exp(-oI1*(t))*np.sin(oR1*(t)+phi1)+\
        A2*np.exp(-oI2*(t))*np.sin(oR2*(t)+phi2)
    return p
#-----------------------------------------------
def findN(a): # number of modes included in the fit 
    n = np.size(a)
    N = 0 if n <4 else 1
    return N 
#--------------------------------------
print('#----------------------------------------------------------')
R      = [float(x) for x in input('# Enter position of observer R : \n').split()]
tshift = [float(x) for x in input('# Enter tshift from tpeak      : \n').split()]
print('#----------------------------------------------------------')
print(f'R      = {R}')
print(f'tshift = {tshift}')
print('#----------------------------------------------------------')
sigma   = -3 
m       = 0 
s       = 5
Nxt     = 1000
potInfo = f'GRaxial{sigma}_k1'##f'GRaxial{sigma}'
icInfo  = 'Gaussian_m{}_s{}'.format(m, s) 
startLoading = time.time()
file    = '../outputs/{}_{}_Nxt{}.npz'.format(potInfo, icInfo, Nxt)

data    = np.load(file)
psi     = data['psi']
t       = data['t']
xt      = data['xt']
#meta    = data['metadata']
#-
tf  = t[-1];
xti = xt[0]; xtf = xt[-2]
endLoading = time.time()
print(f't_loading={endLoading-startLoading}')

#----
print(f'# R before = {R}')
NR       = np.where(xt>=R)[0][0]
R        = xt[NR]  
print(f'# R after = {R}')
h        = psi[:, NR]
logh     = np.log10(np.abs(h))
# tepak 
hmax     = np.max(h) 
maxIndex = np.where(h == hmax)
tpeak    = t[maxIndex]
print(f'tpeak={tpeak}')
# restric t  
tmax     = 95
maxIndex = np.where(t>=tmax)[0][0]
tmin     = tpeak  + tshift
minIndex = np.where(t>=tmin)[0][0]
tfit     =    t[minIndex:maxIndex]
hfit     =    h[minIndex:maxIndex]
loghfit  = logh[minIndex:maxIndex]

#----------fit 
initial_guess = [0.05, 0.2, 0.03, np.pi] 
params, cov = curve_fit(modelQNMorg, tfit, hfit, initial_guess)
N = findN(initial_guess)
print("#A,    omegaR,     omegaI,    phi")
print(params)


#------------------------
#-- logplot
fig=plt.figure(); 
ax =fig.add_subplot(1,1,1)
ax.plot(t, logh)
ax.set_xlabel('t')
ax.set_ylabel ('log|h(t)|')
ax.set_title(f'@ xt = {R}')
ax.axvline(x=tpeak, color="r")
ax.axvline(x=tmin, color="black")
ax.axvline(x=tmax, color="black")
ax.grid()
fig.savefig(f'../figures/fits/h_R{int(R)}.pdf')

#-- |h(t)|
#plot2d.plot2d(t,   h , " ", 't', 'h(t)')


#--- plot to compare
fig=plt.figure()
ax = fig.add_subplot(2,1,1)
ax.plot(tfit, hfit, "b", label='solution')
ax.plot(tfit, modelQNMorg(tfit, *params), '-r', label='fit') 
ax.grid()
ax.set_title(f'N={N}, @ xt= {R}, tshift={tshift}')
ax.set_xlabel('t')
ax.set_ylabel('h(t)')
#--- residual
ax1 = fig.add_subplot(2,1,2)
ax1.plot(tfit, np.abs(hfit-modelQNMorg(tfit, *params)), '--b')
ax1.grid()
ax1.set_title('residual')
ax1.set_xlabel('t')
ax1.set_ylabel(r'$|h_{fit}-h_{model}|$')
ax.legend()
fig.savefig(f'../figures/fits/compareResidual_{int(R)}.pdf')

#--- log plot to compare
fig=plt.figure()
ax = fig.add_subplot(2,1,1)
ax.plot(tfit, loghfit, "b", label='solution')
ax.plot(tfit, np.log10(np.abs(modelQNMorg(tfit, *params))), '-r', label='fit') 
ax.grid()
ax.set_title(f'N={N}, @ xt = {R}, tshift={tshift}')
ax.set_xlabel('t')
ax.set_ylabel('log|h(t)|')
#--- residual
ax1 = fig.add_subplot(2,1,2)
ax1.plot(tfit, np.abs(np.log10(np.abs(hfit-modelQNMorg(tfit, *params)))), '--b')
ax1.grid()
ax1.set_title('residual')
ax1.set_xlabel('t')
ax1.set_ylabel(r'$|log(|h_{fit}|)-log(|h_{model}|)|$')
ax.legend()
fig.savefig(f'../figures/fits/compareResidualLog_{int(R)}.pdf')

#--- copmlex omega plane for different tshift
print(r"#--------- compelex $\omega$ plane ----------------------")
Nof = 10
Dtshift = int((tmax - tpeak[0])//Nof)
tstart  = np.linspace(tpeak, tmax, Nof)

omegaR = []
omegaI = []

# make the plot 
fig, (ax1,ax2) = plt.subplots(2,1)
ax1.set_title(rf'N={N}, Complex plane for $\omega$ for varying tpeak')
ax1.set_xlabel(r'$Re[\omega]$')
ax1.set_ylabel(r'$Im[\omega]$')
ax1.grid(); ax2.grid()
ax2.set_xlabel("tstart")
ax2.set_ylabel(r"$\frac{|\omega_{fit}-\omega_{true}|}{\omega_{true}}$")
ax2.set_yscale("log")
# make the fits 
omega0=0.74927369
for i in range(0, Nof-1):
    # determine fit limits
    minIndex = np.where(t>=tstart[i])[0][0]
    tfit     =    t[minIndex:maxIndex]
    hfit     =    h[minIndex:maxIndex]
    loghfit  = logh[minIndex:maxIndex]
    # fit 
    initial_guess = [0.05, 0.2, 0.03, np.pi] 
    params, cov   = curve_fit(modelQNMorg, tfit, hfit, initial_guess)
    N             = findN(initial_guess)
    omegaR.append(params[1])
    omegaI.append(params[2])
    # print useless info
    print(i)
    print(minIndex)
    print(f"# tstart = {tstart[i]} ")
    print("#A,           omegaR,         omegaI,        phi")
    print(params)
    print("#--------------------------------")
    # fill the plot
    ax1.scatter(omegaR[i], omegaI[i], s=80, marker='o', edgecolors='black', facecolors=('blue',i/(Nof-2)))
    ax2.scatter(tstart[i], np.abs((omegaR[i]-omega0))/omega0,\
                    s=70, marker='o', edgecolors='black', facecolors=('blue', i/(Nof-2)))
plt.tight_layout()        
fig.savefig(f'../figures/fits/complexOmega_R{int(R)}.pdf')
#----
plt.show()