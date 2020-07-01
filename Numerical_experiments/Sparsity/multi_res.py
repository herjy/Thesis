import numpy as np
import pyfits as pf
import subprocess as sp
import os
import time
import shlex
import matplotlib.pyplot as plt

def mr_transform(Imag, MR_File_Name=0, Opt='-t 2', pad0 = 0, combine = 0):
    
    vsize = np.shape(Imag)
    
    if np.size(vsize) ==3:
        mn = vsize[np.int_(np.where(vsize != vsize[1]))]

        if mn == vsize[0]:
                w0 = mr_transform(Imag[0,:,:])
        else:
                w0 = mr_transform(Imag[:,:,0])
        lvl = np.min(np.shape(w0))
        wave = np.zeros([vsize[1], vsize[1],lvl, mn])
        for h in np.linspace(0,mn-1, mn):
            if mn == vsize[0]:
                w=mr_transform(Imag[h,:,:])
                wave[:,:,:,h] = np.transpose(w)#,(vsize[1],vsize[1],lvl))
                
            else:
                wave[:,:,:,h] = mr_transform(Imag[:,:,h])
        return wave
    
    if pad0 != 0:
        pad = np.zeros([pad0,pad0])
        pad[(pad0-vsize[0])/2:(pad0+vsize[0])/2,(pad0-vsize[0])/2:(pad0+vsize[0])/2] = Imag
        Imag = pad
        
    if MR_File_Name !=0:
        filename = MR_File_Name
        connard = 0
    else:
        filename = 'xx_temp'
        connard = 1
    filename = filename + '.mr'
    NameImag='xx_imag.fits'

    hdu = pf.PrimaryHDU(Imag)
   
    pf.HDUList([hdu]).writeto(NameImag, clobber = True)


    com = './mr_transform ' + Opt + '  ' + NameImag + '  ' +  filename
    sp.Popen(com,stdout=sp.PIPE, shell = True).communicate()
    os.remove(NameImag)

    #time.sleep(10)

    if MR_File_Name==0:
        MR_File_Name = './xx_temp.fits'
    sp.call('mv '+filename+' '+MR_File_Name,stdout=sp.PIPE, shell = True)

    Cube = pf.open(MR_File_Name)[0].data
    
    if connard == 1:
        os.remove('./xx_temp.fits')

    
    if combine !=0:
        s0,n1,n2 = np.shape(Cube)
        cucube = np.zeros([n1,n2,(s0-1)/3+1])
        s = s0+0
        
        cucube[:,:,(s0-1)/3] = Cube[s0-1,:,:]
        for k in np.linspace(0,(s0-1)/3-1,(s0-1)/3):
            cucube[:,:,k] = np.sum(Cube[k*3:(3+k*3),:,:],0)/np.max(np.max(np.sum(Cube[k*3:(3+k*3),:,:],0)))
        Cube = cucube
    return Cube

def mr_filter(Imag, MR_File_Name=0, Opt='-f 10 -t 1'):

    vsize = np.shape(Imag)

    if MR_File_Name !=0:
        filename = MR_File_Name
        connard = 0
    else:
        filename = 'xx_temp'
        connard = 1
        
    filename = filename + '.mr'

    NameImag='xx_imag.fits'


    hdu = pf.PrimaryHDU(Imag)
    pf.HDUList([hdu]).writeto(NameImag, clobber = True)

    com = './mr_filter ' + Opt + ' ' + NameImag + ' ' +  filename
    p=sp.Popen(com,stdout=sp.PIPE, shell = True).communicate()

    
    os.remove(NameImag)

    if MR_File_Name==0:
        MR_File_Name = 'xx_temp.fits'
    sp.call('mv '+filename+' '+MR_File_Name,stdout=sp.PIPE, shell = True)

    Cube = pf.open(MR_File_Name)[0].data
    
    if connard == 1:
        os.remove('xx_temp.fits')

    return Cube


def mr_gmca(Imag, MR_File_Name=0, Opt='-t28 -N1 -v -S2'):

    vsize = np.shape(Imag)

    if MR_File_Name !=0:
        filename = MR_File_Name
        connard = 0
    else:
        filename = 'xx_temp.fits'
        connard = 1

    NameImag='xx_imag.fits'

    hdu = pf.PrimaryHDU(Imag)
    pp=pf.HDUList([hdu]).writeto(NameImag, clobber = True)
    
    com = 'mr_gmca' + ' '+ Opt + ' ' + NameImag + ' ' +  filename
    p=sp.call(com,stdout=sp.PIPE, shell = True)#.communicate()

    
    os.remove(NameImag)

    Cube = pf.open(filename)[0].data
    
    if connard == 1:
        os.remove('xx_temp.fits')

    return Cube

    

    
