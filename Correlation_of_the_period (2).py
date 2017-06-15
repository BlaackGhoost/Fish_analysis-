# -*- coding: utf-8 -*-
"""
Created on Wed May  3 13:24:13 2017

@author: Amina
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
from scipy.fftpack import rfft, irfft, fftfreq
from scipy import signal

odata = np.load("Dexter_6.npy")
m1 = []
m2 = []
m3 = []
mp1=[]
s = len(odata) /6
p = 6
#colors = intertools.cycle(["r", "b" , "g" , "y"])

for x in range (1,p+1):
    
    start = int((x-1)*s)
    end = int(x*s)
    test=odata
    test2 = []
    for i in range(start,end):
        test2.append(test[i])
    data = test2
    
    sample_freq = 41666666.666 
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
    thresh=1.2*np.std(EODs)   # 0
    EODs2=np.roll(EODs,-1)  #shifted dummy time series
    cyclei=np.array(np.where((EODs<thresh) & (EODs2> thresh)))[1].T
    EODcross=time[cyclei[1:len(cyclei)-1]]
    EODperiod=np.diff(EODcross,axis=0)
    EODfreq = 1./EODperiod
    aver_period = np.average(EODperiod)
    aver_freq = np.average(EODfreq)
    print('average frequency = ', str(aver_freq))
    m1=np.concatenate((m1,EODfreq),axis=0)
    fs= sample_freq
    t= time
    analytic_signal = hilbert(EODs) # instead of EODs
    amplitude_envelope = np.abs(analytic_signal)
    instantaneous_phase = np.unwrap(np.angle(analytic_signal)) / (2.0*np.pi)
    instant_mod = np.mod(instantaneous_phase, 1)
    instantaneous_frequency = np.diff(instantaneous_phase)  * fs
    aveg_freq_H = np.mean(instantaneous_frequency)
    instant_mod=instant_mod-np.mean(instant_mod)
    instant_mod=instant_mod/np.std(instant_mod)
    thresh2=1.2*np.std(instant_mod)
    instant_mod2=np.roll(instant_mod,-1)
    cyclei3=np.array(np.where((instant_mod<thresh2) & (instant_mod2> thresh2)))[1].T
    instant_mod_cross=time[cyclei3[1:len(cyclei3)-1]]
    instant_mod_period=np.diff(instant_mod_cross,axis=0)
    instant_mod_freq = 1./instant_mod_period
    aver_freq_Hilbert = np.average(instant_mod_freq)
    print('Hilbert average frequency = ', str(aver_freq_Hilbert))
    m2=np.concatenate((m2,instant_mod_freq),axis=0)
    mp1=np.concatenate((mp1,EODperiod),axis=0)

#    colors = "bgrcmykw"
#    color_index= 0
    m11 = (1/m1)*10e2
    m12 = np.roll(m11,1) # shift by one 
    fig3 = plt.figure('correlation')
    plt.clf()
    #cl = plt.cm.coolwarm(m1)
    #plt.plot(m11,m12,'.',c = 'red')
    color = np.sort(m1)
    plt.scatter(m11, m12, c=color, marker='.')
    plt.xlabel( 'time in miliseconds')
    plt.ylabel( 'period shifted by one in miliseconds ')
    x = x+1
# MAIN     

j=0
av1 = []
av2 = []

hist1 = []
hist2 = []


average1 = np.average(m1)
average2 = np.average(m2)

print('Average frequency for all = ', str(average1))
print('Hilbert average frequency for all = ', str(average2))


    


plt.figure('EOD plot')
plt.clf()   
plt.plot(time,EODv[0,:],'.')
plt.plot(time,EODs[0,:], linewidth = 2)
plt.plot([time[0], np.max(time)], [thresh, thresh])
plt.show

fig = plt.figure('fish')
ax0 = fig.add_subplot(211)
ax0.plot(t, EODs[0,:], label='EODs')
ax0.plot(t, amplitude_envelope[0,:], label='envelope')
ax0.set_xlabel("time in seconds")
ax0.legend()
ax1 = fig.add_subplot(212)
ax1.plot(t[1:],instantaneous_frequency[0,:], label = 'frequency' )

ax1.set_xlabel("time in seconds")



binmin= min(m1)
binmax= max(m1)
fig2 = plt.figure('histogram_tresh')
plt.hist(m1, bins=40)#, range=(742,758))
plt.title("histogram_tresh")
plt.xlabel( 'frequency (Hz)')
plt.ylabel( 'no of values')
plt.xlim(binmin,binmax)
plt.show


binmin1= min(instant_mod_freq[1:])
binmax1= max(instant_mod_freq[1:])
fig1 = plt.figure('histogram_mod')
plt.hist(m2, bins = 40)#, range=(742,758))
plt.title("histogram_mod")
plt.xlabel( 'frequency (Hz)')
plt.ylabel( 'no of values')
plt.xlim(binmin1,binmax1)

std_frequency = np.std(m11)
mean_frequency = np.mean(m11)
cv_frequency = std_frequency/mean_frequency 
print('standard deviation frequency  = ', str(std_frequency))
print('mean_frequency = ', str(mean_frequency))
print('CV_frequency = ', str(cv_frequency))

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

signal2 = np.squeeze(mp1)

ml = int(len(signal2)/2)
fig = plt.figure('XCORR')
ax1 = fig.add_subplot(211)
ax1.xcorr(signal2, signal2, usevlines=True, maxlags=ml, normed=True, lw=0.5)
ax1.grid(True)
ax1.axhline(0, color='black', lw=2)

ax2 = fig.add_subplot(212, sharex=ax1)
ax2.acorr(signal2 , usevlines=True, normed=True, maxlags=ml, lw=0.5)
ax2.grid(True)
ax2.axhline(0, color='black', lw=2)

plt.show()


m13 = mp1
m12 = np.roll(m13,1) # shift by one 
fig3 = plt.figure('correlation')
plt.clf()
#cl = plt.cm.coolwarm(m1)
#plt.plot(m11,m12,'.',c = 'red')
color = np.sort(m13)
plt.scatter(m13, m12, c=color, marker='.')
plt.xlabel( 'time in miliseconds')
plt.ylabel( 'period shifted by one in miliseconds ')
x = x+1


f, Pxx_den = signal.welch(EODv, sample_freq, scaling='spectrum', nperseg=1000000)
Pxx_den = np.squeeze(Pxx_den)
plt.clf() 
plt.semilogy(f, Pxx_den)
#plt.ylim([0.5e-3, 1])
plt.xlim(100,10000)
plt.xlabel('frequency [Hz]')
plt.ylabel('Linear spectrum [V RMS]')
plt.show()

plt.figure('frequency with time')
plt.clf() 
t = np.arange(0,(len(m11))*period,period)

plt.plot(t, m11,'.')

#m11 = m1 
#m12 = np.roll(m1,1)
#
#fig3 = plt.figure('negative auto correlation')
#plt.clf()
#plt.plot(m11,m12,'.')