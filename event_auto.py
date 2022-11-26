# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 17:29:27 2022

@author: Administrator
"""

import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pyedflib import highlevel
import numpy as np
import scipy
from scipy.signal import savgol_filter
import xlrd
import pandas
import copy
import heapq
from sampen import sampen2



name = "F:/osa data/EDF/{9766EDA0-2A64-4B61-BC16-1A776133F731}.SLP.edf"
arousal_path = "F:/osa data/Label_Arousal/{9766EDA0-2A64-4B61-BC16-1A776133F731}.SLP_event.xlsx"
# edf
#==============================================================================        
signals, signal_headers, header = highlevel.read_edf(name)        
dic = highlevel.read_edf_header(name)
channels = dic['channels']
print(channels)
#讀取並印出所有channels

nPress = signals[11]  #讀取PSG中的NPress
therm = signals[12] #讀取PSG中的therm
chest = signals[13] #胸腔變化
abdo = signals[14] #腹腔變化
rawSpO2 = signals[15] #血氧濃度變化


fsamp=25 # 資料點/每秒
sleep_lengh= len(nPress) #設定資料總長
time=np.arange(0,(1/fsamp)*sleep_lengh,(1/fsamp)) #與資料點相匹配的時間軸
nPress_sm = scipy.signal.savgol_filter(nPress,81,3)
all_lengh= len(nPress)
fsamp = 25;  # Sample rate

# ==================================================================
# spo2
SpO2 = np.full((all_lengh), 0.2)        
if len(rawSpO2)<all_lengh:
    for i in range(0,len(rawSpO2)):
         SpO2[i*25:i*25+25] = rawSpO2[i] # 擴增血氧濃度資料點
else:
     SpO2=rawSpO2
for i in range(0,all_lengh):
    if SpO2[i]<SpO2[i-1]-30:
        SpO2[i]=SpO2[i-1]
SpO2_base = np.full((all_lengh),np.mean(SpO2[0:fsamp*5*60]))
SpO2_base_va = np.mean(SpO2[0:fsamp*5*60])   
SpO2[all_lengh-(5*fsamp*60):all_lengh] = SpO2_base_va   
SpO2[0:5*fsamp*60] = SpO2_base_va 
time=np.arange(0,(1/fsamp)*all_lengh,(1/fsamp))

# ==================================================================
# arousal import
ar_o2_event = np.full((all_lengh), 0)
arousal = pandas.read_excel(arousal_path)
Astart = arousal['Start'][:]  # arousal開始時間
Adur = arousal['Duration'][:]  # arousal開始時間
for i in range(0,len(Astart)):
    ar_o2_event[round(Astart[i]*25):round((Astart[i]+Adur[i])*25)]=[1]
# ==================================================================
# ventilation        
nPress_sign = np.sign(nPress_sm)
nPress_pos = copy.deepcopy(nPress)
nPress_pos[nPress_pos<0]=0
cp1=[]
cp2=[]
for i in range(0,all_lengh-1):   
    if nPress_sign[i+1]-nPress_sign[i]==2:
        k=i+1
        cp1.append(k)
        for j in range(k,all_lengh-1):
            if nPress_sign[j+1]-nPress_sign[j]==2:
                m=j
                cp2.append(m)
                break
            elif j+1==all_lengh-1:
                m=j+1
                cp2.append(m)
                break     
            
for i in range(0,len(cp1)):            
    nPress_pos[cp1[i]:cp2[i]+1] = [sum(nPress_pos[cp1[i]:cp2[i]+1])/(cp2[i]-cp1[i]+2)]

Vt = nPress_pos
Vt[Vt<np.mean(Vt)*0.1]=0
# ==================================================================
# baseline
Vt_base = np.full((all_lengh), 0.2)
osa_event = np.full((all_lengh), 0)
Vt_base[0:10*fsamp*60-1] = np.mean(Vt[0:10*fsamp*60-1])
        

begain = 1;
k=1;
while k<=all_lengh-5*fsamp*60:
    if Vt[k+1]<Vt_base[k]*0.1 :
        low_beg = k+1
        for j in range(k+1,all_lengh):
                if Vt[j]>Vt_base[k]*0.1:
                   low_end=j-1;
                   break
        if low_end-low_beg>=9.5*fsamp:
               Vt_base[begain:low_end+1] = np.mean(Vt[begain:k+1])
               osa_event[low_beg:low_end+1]=1
               k=j+10*fsamp
               begain=j
               Vt_base[k-10*fsamp:k+1]=np.mean(Vt[k-10*fsamp:k+1])
        else :
            Vt_base[begain:low_end+1] = np.mean(Vt[begain:low_end+1])
            k = low_end+1
    elif Vt[k+1]<=Vt_base[k]*0.7 :
        low_beg = k+1
        for j in range(k+1,all_lengh):
            if Vt[j]>Vt_base[k]*0.4:
               low_end=j-1;
               break
        if low_end-low_beg>=9.5*fsamp:
               Vt_base[begain:low_end+1] = np.mean(Vt[begain:k+1])                       
               for j in range(low_beg,low_end+20*fsamp):
                   if SpO2[i]<SpO2_base[i]*0.97:
                       ar_o2_event[i] = 1                      
               apnea_lengh = 0
               for j in range(low_beg,low_end):
                   if Vt[i]<Vt_base[k]*0.1:                     
                       apnea_lengh=apnea_lengh+1
               if apnea_lengh>=10*fsamp:
                   osa_event[low_beg:low_end+1]=1
                   k=j+10*fsamp
                   begain=j
                   Vt_base[k-10*fsamp:k+1]=np.mean(Vt[k-10*fsamp:k+1])
               elif apnea_lengh<10*fsamp and np.mean(ar_o2_event[low_end:low_end+20*fsamp+1])>0:
                   osa_event[low_beg:low_end+1]=1
                   k=j+10*fsamp
                   begain=j
                   Vt_base[k-10*fsamp:k+1]=np.mean(Vt[k-10*fsamp:k+1])                           
               else :
                   Vt_base[begain:low_end+1] = np.mean(Vt[begain:low_end+1])
                   k = low_end
        else :
            Vt_base[begain:low_end+2] = np.mean(Vt[begain:low_end+2])
            k = low_end+1
    elif Vt[k+1]>Vt_base[k]*0.7 :
        Vt_base[begain:k+fsamp*5]=np.mean(Vt[begain:k+fsamp*5])
        k=k+1
#===================================================================
# event count
event_start=[]
event_end=[]
nu_event = 0
i=0
while i < all_lengh-60*fsamp:                    
    if osa_event[i]==1 and Vt_base[i-1] != 0 :
        event_start.append(i)                
        for j in range(i,all_lengh-1):           
            if osa_event[j]==0:
                nu_event = nu_event +1
                event_end.append(j)
                i=j
                break
    else :
        i=i+1

#============================================================================== 
# figure
ax = plt.subplot(2,1,1)
ax.plot(time/60,nPress,'k')
ax.legend(['npress'])
for i in range(0,nu_event):  
    ax.add_patch(patches.Rectangle((event_start[i]/60/25, -2),(event_end[i]-event_start[i])/60/25,4,
                                           facecolor = 'blue',alpha=0.3,fill=True) ) 
ax2 = plt.subplot(2,1,2, sharex = ax)

ax2.plot(time/60,chest*10,'k')
ax2.plot(time/60,abdo*10,'b')
ax2.legend(['chest','abdo'])