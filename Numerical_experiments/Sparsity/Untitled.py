import pyfits as pf
import matplotlib.pyplot as plt
import numpy as np
import multi_res as mr
import math
import matplotlib.cm as cm
from MuSCADeT import wave_transform as mw
import scipy.ndimage.filters as sc

def lin_upsample(img):
    n1,n2 = np.shape(img)
    upimg = np.zeros((n1*2,n2))
    for i in range(n1):
        upimg[i*2,:] = img[i,:]
    return upimg

def col_upsample(img):
    n1,n2 = np.shape(img)
    upimg = np.zeros((n1,n2*2))
    for i in range(n2):
        upimg[:,i*2] = img[:,i]
    return upimg


def iuwt(W00,w03,w30,w33,w22,w02,w20,w00,w01,w10):
    h = np.array([1.,4.,6.,4.,1.])
    h = h/8.
    ht = h
    i0=1
    i1=0
#    plt.imshow(-W00); plt.show()
    ###premier filtrage
    Wi00 = col_upsample(W00)
    wi03 = col_upsample(w03)
    wi30 = col_upsample(w30)
    wi33 = col_upsample(w33)

    #filtrage colonne 
    Wi00 = sc.convolve1d(Wi00,h,axis = i0, mode = 'constant', cval = 0)
    wii03 = sc.convolve1d(wi03,h,axis = i0, mode = 'constant', cval = 0)
    wi03 = wi03-sc.convolve1d(wii03,h,axis = i0, mode = 'constant', cval = 0)
    wi30 = sc.convolve1d(wi30,h,axis = i0, mode = 'constant', cval = 0)
    wii33 = sc.convolve1d(wi33,h,axis = i0, mode = 'constant', cval = 0)
    wi33 = wi33-sc.convolve1d(wii33,h,axis = i0, mode = 'constant', cval = 0)

    
    W31 = lin_upsample(Wi00+wi03)
    W32 = lin_upsample(wi30+wi33)

    #Filtrage ligne
    W31 = sc.convolve1d(W31,h,axis = i1, mode = 'constant', cval = 0)
    Wi32 = sc.convolve1d(W32,h,axis = i1, mode = 'constant', cval = 0)
    W32 = W32-sc.convolve1d(Wi32,h,axis = i1, mode = 'constant', cval = 0)

    W3 = W32+W31
    
 #   plt.imshow(-W3); plt.show()
    ###Deuxieme filtrage
    Wi3 = col_upsample(W3)
    wi02 = col_upsample(w02)
    wi20 = col_upsample(w20)
    wi22 = col_upsample(w22)

    #filtrage colonne 
    Wi3 = sc.convolve1d(Wi3,h,axis = i0, mode = 'constant', cval = 0)
    wii02 = sc.convolve1d(wi02,h,axis = i0, mode = 'constant', cval = 0)
    wi02 = wi02-sc.convolve1d(wii02,h,axis = i0, mode = 'constant', cval = 0)
    wi20 = sc.convolve1d(wi20,h,axis = i0, mode = 'constant', cval = 0)
    wii22 = sc.convolve1d(wi22,h,axis = i0, mode = 'constant', cval = 0)
    wi22 = wi22-sc.convolve1d(wii22,h,axis = i0, mode = 'constant', cval = 0)
  
    
    W21 = lin_upsample(Wi3+wi02)
    W22 = lin_upsample(wi20+wi22)

    #Filtrage ligne
    W21 = sc.convolve1d(W21,h,axis = i1, mode = 'constant', cval = 0)
    Wi22 = sc.convolve1d(W22,h,axis = i1, mode = 'constant', cval = 0)
    W22 = W22-sc.convolve1d(Wi22,h,axis = i1, mode = 'constant', cval = 0)

    W2 = W22+W21
#    plt.imshow(-W2); plt.show()
    ###troisieme filtrage
    Wi2 = col_upsample(W2)
    wi01 = col_upsample(w01)
    wi10 = col_upsample(w10)
    wi00 = col_upsample(w00)
    
    #filtrage colonne 
    Wi2 = sc.convolve1d(Wi2,h,axis = i0, mode = 'constant', cval = 0)
    wii01 = sc.convolve1d(wi01,h,axis = i0, mode = 'constant', cval = 0)
    wi01 = wi01-sc.convolve1d(wii01,h,axis = i0, mode = 'constant', cval = 0)
    wi10 = sc.convolve1d(wi10,h,axis = i0, mode = 'constant', cval = 0)
    wii00 = sc.convolve1d(wi00,h,axis = i0, mode = 'constant', cval = 0)
    wi00 = wi00-sc.convolve1d(wii00,h,axis = i0, mode = 'constant', cval = 0)

    
    W11 = lin_upsample(Wi2+wi01)
    W12 = lin_upsample(wi10+wi00)

    #Filtrage ligne
    W11 = sc.convolve1d(W11,h,axis = i1, mode = 'constant', cval = 0)
    Wi12 = sc.convolve1d(W12,h,axis = i1, mode = 'constant', cval = 0)
    W12 = Wi12-sc.convolve1d(Wi12,h,axis = i1, mode = 'constant', cval = 0)

    W1 = W11+W12
#    plt.imshow(-W1); plt.show()
    return W1



    
    
    




c = -1


im0 = pf.open('67P.fits')[0].data
ni1,ni2 = np.shape(im0)
im = np.zeros((4096,4096))
im[:ni1,:] = im0[:,:4096]
n1 = 4096; c = 1


im0 = pf.open('minion2.fits')[0].data
im0 = -im0+255
im0 =im0[::-1,:]
z = np.zeros(im0.shape)
x,y = np.where(z==0)
for i in range(len(x)):
    if math.isnan(im0[x[i],y[i]])==True:
        im0[x[i],y[i]] = 0
n1 =2000
ni1,ni2 = np.shape(im0)
print(ni1,ni2)
im = np.zeros((n1,n1))
im[:ni1,:] = im0[:,:]
hdus = pf.PrimaryHDU(im)
lists = pf.HDUList([hdus])
lists.writeto('new_minion.fits', clobber=True)
c = -1

wim = mr.mr_transform(im, Opt = '-t 14')

sim = np.sort(np.reshape(np.abs(im),im.size))
argsim = np.argsort(np.reshape(np.abs(im),im.size))
argsim = argsim[::-1]
sim = sim[::-1]
swim = np.sort(np.reshape(np.abs(wim),wim.size))
argswim = np.argsort(np.reshape(wim,wim.size))
argswim = argswim[::-1]
swim = swim[::-1]

sim[sim.size/100.:] = 0
swim[swim.size/100.:] = 0


tim = sim+0.
twim = swim+0.

x = np.where(np.zeros(n1*n1)==0)
tim[argsim] = sim[x]
tim = np.reshape(tim,(n1,n1))

x = np.where(np.zeros(n1*n1)==0)
twim[argswim] = swim[x]
twim = np.reshape(twim,(n1,n1))

print(np.size(np.where(tim!=0)),np.size(np.where(twim!=0)))
#high frequencies
w00 = twim[n1/2:,n1/2:]
w10 = twim[n1/2:,:n1/2]
w01 = twim[:n1/2,n1/2:]

#second scale
w22 = twim[n1/4:n1/2,n1/4:n1/2]
w20 = twim[n1/4:n1/2,:n1/4]
w02 = twim[:n1/4,n1/4:n1/2]

#first scale
w33 = twim[n1/8.:n1/4.,n1/8.:n1/4.]
w30 = twim[n1/8.:n1/4.,:n1/8.]
w03 = twim[:n1/8.,n1/8.:n1/4.]
#coarse
W00 = twim[:n1/8.,:n1/8.] 
print(W00.shape,w33.shape,w22.shape,w00.shape)

wimrec = iuwt(W00,w03,w30,w33,w22,w02,w20,w00,w01,w10)


hdus = pf.PrimaryHDU(wimrec)
lists = pf.HDUList([hdus])
lists.writeto('wave_minion.fits', clobber=True)


plt.figure(1)
plt.imshow(im, cmap = 'hsv', interpolation = None)
#plt.title('Original image')
plt.axis('off')
plt.figure(2)
plt.imshow((np.abs(wim))**0.5, cmap = 'hsv', interpolation = None)
#plt.title('Wavelet decomposition')
plt.axis('off')
plt.figure(4)
plt.imshow(np.abs(tim), vmin = np.min(tim[tim>0]), cmap = 'hsv', interpolation = None)
#plt.title('1% brightest pixels')
plt.axis('off')
plt.figure(5)
plt.imshow(np.abs(twim)**0.5, cmap = 'hsv', interpolation = None)
#plt.title('1% highest wavelet coefficients')
plt.axis('off')
plt.figure(6)
plt.imshow(np.abs(wimrec), cmap = 'hsv', interpolation = None)
#plt.title('reconstruction from 1% wavelet coefficients')
plt.axis('off')
plt.show()
