import Haar_deconv as hd
import pyfits as pf
import matplotlib.pyplot as plt
import gaussian as gs
import numpy as np
import pywt
import SLIT



Image = pf.open('0416_Par.fits')[0].data#[-1,:,:]#[2000:-2000, 2000:-2000]
print(Image.shape)
Show = np.copy(Image.T)*0
for i in range(3):
    Show[:, :, i] = Image[i, :, :]

Show0 = Show/np.max(Image+0.02)-np.min(Image/np.max(Image+0.02))
plt.imshow(Show0)
plt.axis('off')
plt.show()

WS = SLIT.tools.wave_transform(Image, 6)
for i in range(WS.shape[0]):
    print(WS[i,:,:].shape)
    plt.figure(i)
    plt.imshow(WS[i,:,:]/np.max(WS[i,:,:]+0.02)-np.min(WS[i,:,:]/np.max(WS[i,:,:]+0.02)))#, cmap = 'gray', interpolation=None)
    plt.axis('off')
plt.show()


Im_line = np.sort(np.abs(Image).flatten())[::-1]
Wave_line = np.sort(np.abs(WS).flatten())[::-1]


plt.plot(np.array(range(Im_line.size))*100./Im_line.size, Im_line/np.max(Im_line), 'r')
plt.plot(np.array(range(Wave_line.size))*100./Wave_line.size, Wave_line/np.max(Wave_line), 'b')
plt.show()

Image[np.abs(Image)<Im_line[0.1*Im_line.size]] = 0
WS[np.abs(WS)<Wave_line[0.1*Wave_line.size]] = 0

print(WS.shape)


ShowW = np.copy(Image.T)*0
for i in range(3):
    ShowW[:, :, i] = SLIT.tools.iuwt(WS[:, :, :,i])

Show = np.copy(Image.T)*0
for i in range(3):
    Show[:, :, i] = Image[i, :, :]
plt.figure(0)
plt.imshow((Show/np.max(Image+0.02)-np.min(Image/np.max(Image+0.02))))
plt.axis('off')
plt.figure(1)
plt.imshow((ShowW/np.max(Image+0.02)-np.min(Image/np.max(Image+0.02))))
plt.axis('off')
plt.show()



###############################################################

#Image = pf.open('67P.fits')[0].data
#[:3498,:3498]
#n1,n2 = np.shape(Image)
Image = pf.open('0416_Par.fits')[0].data[0,:,:]#Image[:3499,400:3899]*1.
print(Image.shape)
n1,n2 = np.shape(Image)

print(n1,n2)
W = pywt.wavedec2(Image, 'Haar', level =2)
W = hd.Image_Haar(W, n1,n2, visu = 1)


plt.plot(np.array(range(Image.size))*100./Image.size, np.sort(np.abs(Image.flatten()))[::-1]/np.max(Image), 'r')
plt.plot(np.array(range(W.size))*100./W.size, np.sort(np.abs(np.array(W).flatten()))[::-1]/np.max(W), 'b')
plt.show()

plt.figure(0)
plt.imshow((hd.Image_Haar(W, n1,n2, visu = 1)), cmap = 'gray', interpolation = None, vmin = 0)
plt.axis('off')
plt.figure(1)
plt.imshow(Image, cmap = 'gray', interpolation = None, vmin = 0)
plt.axis('off')
plt.show()

