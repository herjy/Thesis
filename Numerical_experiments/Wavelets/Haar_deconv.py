import pyfits as pf
import numpy as np
from scipy import signal as scp
import scipy.ndimage.filters as sc
import matplotlib.pyplot as plt
import pywt
import SLIT


def levels(n1,n2,lvl, k):
    dirac = np.zeros((n1,n2))
    dirac[n1/2,n2/2] = 1
    wave = pywt.wavedec2(dirac, 'Haar', level =lvl)

    level = [wave[0]*0]
    for i in range(lvl):

        l0 = np.sqrt(np.sum(wave[i+1][0]**2))*k+wave[i+1][0]*0
        l1 = np.sqrt(np.sum(wave[i+1][1]**2))*k+wave[i+1][1]*0
        l2 = np.sqrt(np.sum(wave[i+1][2]**2))*k+wave[i+1][2]*0
        level.append((l0,l1,l2))

    return level

def Image_Haar(W, n1,n2, visu = 0):
    X = np.zeros((n1,n2))
    lvl = np.shape(W)[0]
    if visu == 1:
        X[:n1/2**(lvl-1)+1,:n2/2**(lvl-1)+1] = W[0]/np.max(W[0])
    else:
        X[:n1 / 2 ** (lvl - 1) + 1, :n2 / 2 ** (lvl - 1) + 1] = W[0]
    for i in range(lvl-1):
        j = i+1


        if visu == 1:
            X[n1/2**(lvl-j):n1/2**(lvl-j)+np.shape(W[j][2])[0],n2/2**(lvl-j):n2/2**(lvl-j)+np.shape(W[j][2])[1]] = W[j][2]/np.max(np.abs(W[j][2]))
            X[n1 / 2 ** (lvl - j):n1 / 2 ** (lvl - j)+np.shape(W[j][2])[0], : np.shape(W[j][2])[0]] = W[j][1]/np.max(np.abs(W[j][1]))
            X[:np.shape(W[j][2])[0], n1 / 2 ** (lvl - j):n1 / 2 ** (lvl - j)+np.shape(W[j][2])[0]] = W[j][0]/np.max(np.abs(W[j][0]))
        else:
            X[n1/2**(lvl-j):n1/2**(lvl-(j+1)),n2/2**(lvl-j):n2/2**(lvl-(j+1))] = W[j][2]
            X[n1 / 2 ** (lvl - j):n1 / 2 ** (lvl - (j + 1)), : n1 / 2 ** (lvl - j)+k] = W[j][1]
            X[:n1 / 2 ** (lvl - j)+k, n1 / 2 ** (lvl - j):n1 / 2 ** (lvl - (j + 1))] = W[j][0]

    return X

def wave_line(X, lvl):
    n1,n2 = X.shape
    W = [X[:n1/2**(lvl)+1,:n2/2**(lvl)+1]]
    for i in range(lvl):
        j = i
        if i==2:
            k=0
        else:
            k=1
        w0 = X[:n1 / 2 ** (lvl - j)+k, n1 / 2 ** (lvl - j):n1 / 2 ** (lvl - (j + 1))]
        w1 = X[n1 / 2 ** (lvl - j):n1 / 2 ** (lvl - (j + 1)), : n1 / 2 ** (lvl - j)+k]
        w2 = X[n1/2**(lvl-j):n1/2**(lvl-(j+1)),n2/2**(lvl-j):n2/2**(lvl-(j+1))]
        W.append((w0,w1,w2))
    return W

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


def Soft(X, levels, sigma):
    #   Soft thresholding operator
    lvl = np.size(X)
    Xnew = [X[0]]
    for i in range(lvl-1):

        X0 = np.sign(X[i+1][0])*(np.abs(X[i+1][0])-levels[i+1][0]*sigma)
        X0[(np.abs(X[i+1][0])-levels[i+1][0]*sigma)<0] = 0

        X1 = np.sign(X[i+1][1]) * (np.abs(X[i+1][1]) - levels[i+1][1] * sigma)
        X1[(np.abs(X[i+1][1]) - levels[i+1][1] * sigma) < 0] = 0

        X2 = np.sign(X[i+1][2]) * (np.abs(X[i+1][2]) - levels[i+1][2] * sigma)
        X2[(np.abs(X[i+1][2]) - levels[i+1][2] * sigma) < 0] = 0

        Xnew.append((X0,X1,X2))

    return Xnew

def grad(alpha,res, mu, lvl):
    alphanew = np.copy(alpha)*0
    alphanew = [alpha[0] + mu*res[0]]
    for j in range(lvl):
        a0 = alpha[j+1][0] + mu*res[j+1][0]
        a1 = alpha[j + 1][1] + mu * res[j + 1][1]
        a2 = alpha[j + 1][2] + mu * res[j + 1][2]
        alphanew.append((a0,a1,a2))
    return alphanew

def convolve(X,Y):
    return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X*Y))))
def Dconvolve(X,Y):
    return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X/ Y))))

def deconv(Y, H, niter, k, lvl):
    HT = np.rot90(H)
    n1,n2 = np.shape(Y)

    def PSF_apply(i):
        return scp.fftconvolve(i, H, mode='same')

    def PSFT_apply(ii):
        return scp.fftconvolve(ii, HT, mode='same')
    def Phi(i):
        return pywt.wavedec2(i,'Haar', level = lvl)
    def PhiT(i):
        return pywt.waverec2(i, 'Haar')

    sigma = SLIT.tools.MAD(Y)
    level = levels(n1, n2, lvl, k*np.sqrt(np.sum(HT**2)))

    mu =0.0001
    alpha = Phi(Y*0)

    for i in range(niter):
        print(i)
        Halpha = PhiT(alpha)

        Res = Phi(PSFT_apply(Y-PSF_apply(Halpha)))

        alpha = grad(alpha,Res, mu, lvl)

        alpha = Soft(alpha, level, sigma)


    return PhiT(alpha)