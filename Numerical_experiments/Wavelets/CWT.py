import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as scp


t = np.linspace(0,20,10000)
t_52 = np.where((5<=t)*(t<10))
f_t = 3*np.cos(2*t)*np.sin(5*t)
f_t[t_52] = np.sin(15*t[t_52])+np.cos(23*t[t_52])
f_t[7000:] += 8
f_t[:3000] = 12*np.exp(-0.4*(t[:3000]-2)**2)-2
f_t[5000:6000] = 0
f_t+=np.sin(0.005*t)

f_t/=np.max(f_t)
f_t-=np.mean(f_t)


def Starlet(t):
    return 1./12.*(np.abs(t-2)**3-4*np.abs(t-1)**3+6*np.abs(t)**3-4*np.abs(t+1)**3+np.abs(t+2)**3)

x = np.linspace(-6,6,1000)
plt.figure(0)
plt.plot(x,Starlet(x), 'k', linewidth = 5)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.figure(1)
plt.plot(x,Starlet(x)-1./4.*Starlet(x/2), 'k', linewidth = 5)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.show()

def Mexican(a,tau,t):
    return (1-((t-tau)/a)**2)*np.exp(-(t-tau)**2/(2*a**2))

def Haar(a,tau,t):
    if (t>=(tau-0.5*a))*(t<(tau)):
        return 1./a
    elif (t>=(tau))*(t<(tau+0.5*a)):
        return -1./a
    else:
        return 0

def Morlet(a,tau,t):
    return 1./np.sqrt(2*np.pi)*np.exp(-(t-tau)**2/(2*a**2))*np.cos(2*np.pi*(t-tau)/a)

def trace(f,t):
    F = [f(it) for it in t]
   # plt.plot(F); plt.show()
    return F

def transform_a(Phi,t,s,a):
    Hat_s = s#np.fft.fft(s)
    def Wave(t):
        return Phi(a,10,t)

    Daughter = (trace(Wave,t))
  #  plt.plot(Daughter); plt.show()
    conv = scp.fftconvolve(Hat_s,Daughter, mode = 'same')
    return conv/np.max(conv)##np.fft.ifft(np.fft.fft(Hat_s)*np.conj(np.fft.fft(Daughter)))#

w = 5
a = 2.
def Wave(t):
    return Haar(5., 10, t)
plt.figure(8)
plt.plot(np.array(trace(Wave,t))*4,'k', linewidth = w)
plt.title('Haar\'s wavelet', fontsize = 45)
plt.ylim(-1.1,1.1)
plt.axis('off')
def Wave(t):
    return Morlet(np.float(a), 10, t)
plt.figure(9)
plt.plot(np.array(trace(Wave,t))*2,'k', linewidth = w)
plt.title('Morlet\' wavelet', fontsize = 45)
plt.ylim(-1.1,1.1)
plt.axis('off')
def Wave(t):
    return Mexican(np.float(a), 10, t)
plt.figure(10)
plt.plot(trace(Wave,t),'k', linewidth = w)
plt.title('Mexican hat wavelet', fontsize = 45)
plt.ylim(-1.1,1.1)
plt.axis('off')
plt.show()

size = 1000
wavelet_Haar = np.zeros((size, t.size))
count = 0
for a in np.linspace(2./1000.,2.5,size):
    wavelet_Haar[count, :] = np.real(transform_a(Haar, t, f_t, a))
    count+=1

wavelet_Morlet = np.zeros((size, t.size))
count = 0
for a in np.linspace(2./1000.,2.5,size):
    wavelet_Morlet[count, :] = np.real(transform_a(Morlet, t, f_t, a))
    count+=1

wavelet_Mexican = np.zeros((size, t.size))
count = 0
for a in np.linspace(2./1000.,2.5,size):
    wavelet_Mexican[count, :] = np.real(transform_a(Mexican, t, f_t, a))
    count+=1


f = 20
plt.figure(1)
plt.title('Signal', fontsize = 25)
plt.plot(t, f_t,'k', linewidth = w)
plt.axis('off')
plt.figure(2)
plt.imshow(wavelet_Haar, cmap = 'gray')
plt.title('Haar scalogram', fontsize = 25)
plt.axis('off')
plt.figure(3)
plt.imshow(wavelet_Morlet, cmap = 'gray')
plt.title('Morlet scalogram', fontsize = 25)
plt.axis('off')
plt.figure(4)
plt.imshow(wavelet_Mexican, cmap = 'gray')
plt.title('Mexican hat scalogram', fontsize = 25)
plt.axis('off')
plt.show()

plt.figure(5)
plt.title('scale 0.001', fontsize = 35)
plt.plot(t, f_t, 'k', label = 'signal', linewidth = w)
plt.plot(t,wavelet_Morlet[1], label = 'Morlet', linewidth = w)
plt.plot(t,wavelet_Mexican[1], label = 'Mexican', linewidth = w)
plt.plot(t,wavelet_Haar[1], label = 'Haar', linewidth = w)
plt.legend(fontsize = f)
plt.axis([0,20,-2,2])
a = plt.axes([0.5, 0.65, 0.32, 0.3])
plt.plot(t, f_t, 'k',  linewidth = w)
plt.plot(t,wavelet_Morlet[1], linewidth = w)
plt.plot(t,wavelet_Mexican[1], linewidth = w)
plt.plot(t,wavelet_Haar[1], linewidth = w)
plt.xlim(13.9,14.1)
plt.ylim(-1.8,1.8)
plt.axis('off')

plt.figure(6)
plt.title('Scale 0.07', fontsize = 35)
plt.plot(t, f_t,'k', label = 'signal', linewidth = w)
plt.plot(t,wavelet_Morlet[20], label = 'Morlet', linewidth = w)
plt.plot(t,wavelet_Mexican[20], label = 'Mexican', linewidth = w)
plt.plot(t,wavelet_Haar[20], label = 'Haar', linewidth = w)
plt.legend(fontsize = f)
plt.axis([0,20,-2,2])
a = plt.axes([0.5, 0.65, 0.32, 0.3])
plt.plot(t, f_t,'k', label = 'signal', linewidth = w)
plt.plot(t,wavelet_Morlet[70], label = 'Morlet', linewidth = w)
plt.plot(t,wavelet_Mexican[70], label = 'Mexican', linewidth = w)
plt.plot(t,wavelet_Haar[70], label = 'Haar', linewidth = w)
plt.xlim(13.8,14.2)
plt.ylim(-1.8,1.8)
plt.axis('off')

plt.figure(7)
plt.title('Scale  0.1', fontsize = 35)
plt.plot(t, f_t, 'k', label = 'signal', linewidth = w)
plt.plot(t,wavelet_Morlet[-1], label = 'Morlet', linewidth = w)
plt.plot(t,wavelet_Mexican[-1], label = 'Mexican', linewidth = w)
plt.plot(t,wavelet_Haar[-1], label = 'Haar', linewidth = w)
plt.legend(fontsize = f)
plt.axis([0,20,-2,2])
plt.show()


