

# -*- coding: utf-8 -*-
"""
Created on Fri May 12 10:15:33 2017

@author: Default
"""

import numpy as np
import pandas as pd
from scipy.signal import gaussian
from scipy.ndimage import filters
import matplotlib.pyplot as plt
from math import sqrt
from scipy.signal import hilbert
import os 
from mpl_toolkits.mplot3d import Axes3D
import random 
from scipy.signal import spectrogram


odata = np.fromfile("Lisa609", dtype='i2').reshape(-1, 2).T

odata = odata[0,:]


m1 = []
m11=[]
mp=[]
mp1=[]
s = len(odata) /60
p = 11
l = p-4

for x in range(5,p+1):
    
    start = int((x-1)*s)
    end = int(x*s)
    test=odata
    test2 = []
    for i in range(start,end):
        test2.append(test[i])
    data = test2
    
    sample_freq = 10416666.666 
    period = 1.0/sample_freq
    
    print ('sampling period = ', str(period))
    # file ONE : TRESH
    time=np.arange(0,len(data)*period,period)
    EODv=np.asmatrix(data)
    
    EODv = np.array(EODv,dtype='float64')
    EODv=EODv-np.mean(EODv)
    EODv=EODv/np.std(EODv)
    ff=gaussian(1000,300) #(1000,300)
    EODs = filters.convolve1d(EODv[:], ff/ff.sum()) #filters.convolve1d(EODv[0:100000:10], ff/ff.sum())
    thresh=-1.2*np.std(EODs)   # 0
    EODs2=np.roll(EODs,-1)  #shifted dummy time series
    cyclei=np.array(np.where((EODs<thresh) & (EODs2> thresh)))[1].T
    EODcross=time[cyclei[1:len(cyclei)-1]]
    EODperiod=np.diff(EODcross,axis=0)
    EODfreq = 1./EODperiod
    aver_period = np.average(EODperiod)
    aver_freq = np.average(EODfreq)
    print('average frequency = ', str(aver_freq))
    mp.append(EODperiod)
    m1.append(EODfreq)
    mp1=np.concatenate((mp1,EODperiod),axis=0)
    m11=np.concatenate((m11,EODfreq),axis=0)
    print(x)
    
    x = x+1
# MAIN     

mp2 = np.roll(mp1,1) # shift by one 
fig3 = plt.figure('correlation')
plt.clf()
color = np.sort(mp1)
plt.scatter(mp1, mp2, c=color, marker='.')
plt.xlabel( 'period ms')
plt.ylabel( 'period shifted by one in miliseconds ')


j=0
av1 = []
hist1 = []


minl= []
maxl= []
for k in range (l):
    minl.append(min(m1[k]))
    k = k+1
minbin = min(minl)
for k2 in range (l):
    maxl.append(max(m1[k2]))
    k2 = k2+1
maxbin = max(maxl)
for j in range (l):
    average1 = np.average(m1[j])
   
    
    histogram1= np.histogram(m1[j], bins = 20)# range = (minbin,maxbin))
   


    av1.append(average1)
   
    hist1.append(histogram1)


    j = j+1
print(av1)

plt.figure('EOD plot')
plt.clf()   
plt.plot(time,EODv[0,:],'.')
plt.plot(time,EODs[0,:], linewidth = 2)
plt.plot([time[0], np.max(time)], [thresh, thresh])
plt.show

binmin= min(m11)
binmax= max(m11)
fig2 = plt.figure('histogram')
plt.clf()  
plt.hist(m11, bins=20)#, range=(742,758))
plt.title("histogram")
plt.xlabel( 'frequency [Hz]')
plt.ylabel( 'no of values')
plt.xlim(binmin,binmax)
plt.show

std_frequency = np.std(m11)
mean_frequency = np.mean(m11)
cv_frequency = std_frequency/mean_frequency 
print('standard deviation frequency  = ', str(std_frequency))
print('mean_frequency = ', str(mean_frequency))
print('CV_frequency = ', str(cv_frequency))

std_period = np.std(mp1)
mean_period = np.mean(mp1)
cv_period = std_period/mean_period
print('standard deviation period = ', str(std_period))
print('mean_period = ', str(mean_period))
print('CV_period = ', str(cv_period))


#histogram 
#r = np.ones(2)
#fig = plt.figure()
#ax1 = fig.add_subplot(111, projection = '3d')
#ax1.plot(histogram1)

fig2 = plt.figure()
ax1 = fig2.add_subplot(111, projection='3d')
for n in range (l): 
    

    xn = hist1[n]
    x1= xn[1]
    x1 = np.delete(x1,20,0)
    
    xpos = x1
    ypos = [n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n]
    num_elements = len(xpos)
    zpos = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    dx = np.ones(20)
    #dy = np.ones(10)
    dy = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
    dz = histogram1[0]
    ax1.bar3d(xpos, ypos, zpos, dx, dy, dz)
    
    #ax1.set_xlim3d(500,1000)
    n = n+1

ax1.set_ylim3d(0,l) # (0,l) for number of files 
#ax1.set_xlim3d(800,1000)  




plt.figure('spectr')
plt.clf()  
nperseg = 1000000#10e6
f,t ,Sxx = spectrogram (EODv, sample_freq, nperseg= 1000000 , noverlap= nperseg // 4) 
Sxx = np.squeeze(Sxx)
plt.pcolormesh(t,f, Sxx)
#plt.plot(t,f, Sxx)
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.ylim([0, 5000])

plt.show()
    


plt.figure('frequency with time')
plt.clf() 
t = np.arange(0,len(m11)*period,period)

plt.plot(t, m11,'.')

plt.figure('period with time')
plt.clf() 
t = np.arange(0,len(mp1)*period,period)

plt.plot(t, mp1,'.')




