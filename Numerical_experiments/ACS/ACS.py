import numpy as np
import matplotlib.pyplot as plt

F814 = np.loadtxt('wfc_F814W.dat')
F606 = np.loadtxt('wfc_F606W.dat')
F435 = np.loadtxt('wfc_F435W.dat')

spiral1 = np.loadtxt('Gsbc.spec.txt')
spiral1[:,1]/=np.max(spiral1[:,1])
spiral2 = np.loadtxt('Gscd.spec.txt')
spiral2[:,1]/=np.max(spiral2[:,1])
elliptical = np.loadtxt('Geso.spec.txt')
elliptical[:,1]/=np.max(elliptical[:,1])

plt.plot(F814[:,0],F814[:,1], 'm', linewidth = 5, label = 'F814W')
plt.plot(F606[:,0],F606[:,1], 'g', linewidth = 5, label = 'F606W')
plt.plot(F435[:,0],F435[:,1], 'c', linewidth = 5, label = 'F435W')
plt.plot(spiral1[:,0]*(1+2.381),spiral1[:,1], 'b', linewidth = 5, label = 'Spirals')
plt.plot(spiral2[:,0]*(1+2.381),spiral2[:,1], 'b', linewidth = 5)
plt.plot(elliptical[:,0]*(1+0.44),elliptical[:,1], 'r', linewidth = 5, label = 'Elliptical')
plt.xticks(size = 30)
plt.yticks(size = 30)
plt.axis([3000,11000, 0,1.1])
plt.xlabel('Wavelength in Angstroms', fontsize = 35)
plt.ylabel('Dimensionless thoughput', fontsize = 35)
plt.legend(fontsize = 20)
plt.show()