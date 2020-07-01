import numpy as np
import matplotlib.pyplot as plt


t = np.linspace(0,1000, 10000)

a = 0.4
b = 8
c = a*np.sqrt(1.+b**2/a**2)
theta = -0 *np.pi/180
x = a*(t**2+1)/(2*t)
y = b*(t**2-1)/(2*t)

X = (x)*np.cos(theta) - y*np.sin(theta)
Y = (x)*np.sin(theta) + y*np.cos(theta)

XSun = 0#(c)*np.cos(theta)
YSun = -c/3.#(c)*np.sin(theta)

sel = (Y>-45)*(Y<40)
Xsel = Y[sel]
Ysel = -X[sel]-20

Xstar = np.min(Xsel)
Ystar = Ysel[Xsel==np.min(Xsel)]
print(Xstar, Ystar)
Xobs = np.max(Xsel)
Yobs = Ysel[Xsel==np.max(Xsel)]

plt.plot(Xsel,Ysel,'r', linewidth = 8)
#plt.plot(0,-25, 'oy', markersize = 108, markeredgecolor='y')
#plt.plot(Xstar, Ystar, 'oy', markersize = 108, markeredgecolor='y')
#plt.plot(Xobs,Yobs, 'oy', markersize = 108, markeredgecolor='y')
#plt.xlim([-10,10])
#plt.ylim([-10,10])
plt.axis('off')
plt.show()