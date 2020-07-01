import numpy as np
import matplotlib.pyplot as plt


t = np.linspace(0,2*np.pi,100)

betax = np.sin(t)
betay = np.cos(t)


def lensing(kappa, gamma1, gamma2):
    thetax = 1./((1.-kappa)**2-gamma1**2-gamma2**2)*((1-kappa+gamma1)*betax +gamma2*betay)
    thetay = 1./((1.-kappa)**2-gamma1**2-gamma2**2)*((1-kappa-gamma1)*betay +gamma2*betax)

    return thetax, thetay



plt.plot(betax, betay, 'k', linewidth = 4, label = 'Unit circle')
thetax, thetay = lensing(0.4,0.,0.)
plt.plot(thetax,thetay, 'y', linewidth = 4, label = '$\kappa = 0.4, \gamma_1 = 0, \gamma_2 = 0$')
thetax, thetay = lensing(0.4,1,0.)
plt.plot(thetax,thetay, 'c', linewidth = 4, label = '$\kappa = 0.4, \gamma_1 = 1, \gamma_2 = 0$')
thetax, thetay = lensing(0.4,-1,0.)
plt.plot(thetax,thetay, 'b', linewidth = 4, label = '$\kappa = 0.4, \gamma_1 = -1, \gamma_2 = 0$')
thetax, thetay = lensing(0.4,0.,1)
plt.plot(thetax,thetay, 'r', linewidth = 4, label = '$\kappa = 0.4, \gamma_1 = 0, \gamma_2 = 1$')
thetax, thetay = lensing(0.4,0.,-1)
plt.plot(thetax,thetay, 'm', linewidth = 4, label = '$\kappa = 0.4, \gamma_1 = 0, \gamma_2 = -1$')
plt.xlabel('$\\theta_x$',fontsize = 35)
plt.ylabel('$\\theta_y$',fontsize = 35)
plt.axis([4,-4,4,-4])
plt.legend(fontsize = 20)
plt.show()