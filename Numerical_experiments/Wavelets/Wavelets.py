#coding: utf-8
import numpy as np
import matplotlib.pyplot as plt





tau = np.linspace(-5,5,1000).reshape(1,1000)

def Gabor(nu,tau):
    gauss = np.exp(-np.pi * (-tau) ** 2)
    Real = gauss*np.cos(np.dot(nu.T,tau))
    Imag = gauss*np.sin(np.dot(nu.T,tau))
    return Real, Imag

nu = np.linspace(np.pi/2.,4*np.pi,4).reshape(1,4)
R,I = Gabor(nu, tau)

for n in range(nu.size):
    plt.figure(0)
    plt.plot(tau.T, R[n,:], linewidth = 8, label = '$\\nu$ ='+str( np.round(nu[:,n][0]/np.pi,decimals = 1))+'$\pi$')
    plt.axis([-2,2,-1.1,1.1])
    plt.legend(fontsize = 25)
    plt.axis('off')
    plt.figure(1)
    plt.plot(tau.T, I[n,:], linewidth = 8, label = '$\\nu$ ='+str(np.round(nu[:,n][0]/np.pi,decimals = 1))+'$\pi$')
    plt.axis([-2,2,-1.1,1.1])
    plt.axis('off')
    plt.legend(fontsize = 25)
plt.show()
