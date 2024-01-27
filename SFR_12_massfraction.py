from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from scipy.stats import binned_statistic
from scipy import stats
import numpy as np           # to define our table
import statistics
import matplotlib.mlab
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math                  
import sys

c1=7.301
c5=7.756
c6=8.307
c7=8.757
c8=9.006
c9=9.477
c10=9.760
c11=10
C12=10.146



F1M1B_li=[]
F1M1NB_li=[]
F1M2B_li=[]
F1M2NB_li=[]
F1M3B_li=[]
F1M3NB_li=[]
F1M4B_li=[]
F1M4NB_li=[]


F5M1B_li=[]
F5M1NB_li=[]
F5M2B_li=[]
F5M2NB_li=[]
F5M3B_li=[]
F5M3NB_li=[]
F5M4B_li=[]
F5M4NB_li=[]

F6M1B_li=[]
F6M1NB_li=[]
F6M2B_li=[]
F6M2NB_li=[]
F6M3B_li=[]
F6M3NB_li=[]
F6M4B_li=[]
F6M4NB_li=[]

F7M1B_li=[]
F7M1NB_li=[]
F7M2B_li=[]
F7M2NB_li=[]
F7M3B_li=[]
F7M3NB_li=[]
F7M4B_li=[]
F7M4NB_li=[]

F8M1B_li=[]
F8M1NB_li=[]
F8M2B_li=[]
F8M2NB_li=[]
F8M3B_li=[]
F8M3NB_li=[]
F8M4B_li=[]
F8M4NB_li=[]

F9M1B_li=[]
F9M1NB_li=[]
F9M2B_li=[]
F9M2NB_li=[]
F9M3B_li=[]
F9M3NB_li=[]
F9M4B_li=[]
F9M4NB_li=[]

F10M1B_li=[]
F10M1NB_li=[]
F10M2B_li=[]
F10M2NB_li=[]
F10M3B_li=[]
F10M3NB_li=[]
F10M4B_li=[]
F10M4NB_li=[]

F11M1B_li=[]
F11M1NB_li=[]
F11M2B_li=[]
F11M2NB_li=[]
F11M3B_li=[]
F11M3NB_li=[]
F11M4B_li=[]
F11M4NB_li=[]

F12M1B_li=[]
F12M1NB_li=[]
F12M2B_li=[]
F12M2NB_li=[]
F12M3B_li=[]
F12M3NB_li=[]
F12M4B_li=[]
F12M4NB_li=[]


d='disk_sfr_12_massfraction.dat'
data=np.loadtxt(d,dtype={'names':('WINGS', 'Weights','SFR1', 'sfr_5','sfr_6', 'sfr_7','sfr_8', 'sfr_9','sfr_10', 'sfr_11','sfr_12','M1','M5','M6','M7','M8','M9','M10','M11','M12','Mtotal','F1','F5','F6','F7', 'F8','F9','F10','F11', 'F12','Bar','logM' ),'formats': ('O','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f')} ,skiprows=1)


SFR1=data['SFR1']
sfr5=data['sfr_5']
sfr6=data['sfr_6']
sfr7=data['sfr_7']
sfr8=data['sfr_8']
sfr9=data['sfr_9']
sfr10=data['sfr_10']
sfr11=data['sfr_11']
sfr12=data['sfr_12']
M1=data['M1']
M5=data['M5']
M6=data['M6']
M7=data['M7']
M8=data['M8']
M9=data['M9']
M10=data['M10']
M11=data['M11']
M12=data['M12']
Mtotal=data['Mtotal']
F1=data['F1']
F5=data['F5']
F6=data['F6']
F7=data['F7']
F8=data['F8']
F9=data['F9']
F10=data['F10']
F11=data['F11']
F12=data['F12']
Bar=data['Bar']
Mass=data['logM']




for i in range(len(Mass)):
    if 9.5 < Mass[i] < 10.0 and Bar[i] == 1:
        F1M1B=F1[i]
        F5M1B=F5[i]
        F6M1B=F6[i]
        F7M1B=F7[i]
        F8M1B=F8[i]
        F9M1B=F9[i]
        F10M1B=F10[i]
        F11M1B=F11[i]
        F12M1B=F12[i]
        F1M1B_li.append(F1M1B)
        F5M1B_li.append(F5M1B)
        F6M1B_li.append(F6M1B)
        F7M1B_li.append(F7M1B)
        F8M1B_li.append(F8M1B)
        F9M1B_li.append(F9M1B)
        F10M1B_li.append(F10M1B)
        F11M1B_li.append(F11M1B)
        F12M1B_li.append(F12M1B)

        
        
        
        
    if 9.5 < Mass[i] < 10.0 and Bar[i] == 0:
        F1M1NB=F1[i]
        F5M1NB=F5[i]
        F6M1NB=F6[i]
        F7M1NB=F7[i]
        F8M1NB=F8[i]
        F9M1NB=F9[i]
        F10M1NB=F10[i]
        F11M1NB=F11[i]
        F12M1NB=F12[i]
        F1M1NB_li.append(F1M1NB)
        F5M1NB_li.append(F5M1NB)
        F6M1NB_li.append(F6M1NB)
        F7M1NB_li.append(F7M1NB)
        F8M1NB_li.append(F8M1NB)
        F9M1NB_li.append(F9M1NB)
        F10M1NB_li.append(F10M1NB)
        F11M1NB_li.append(F11M1NB)
        F12M1NB_li.append(F12M1NB)

        
        
    if 10.0 < Mass[i] < 10.5 and Bar[i] == 1:
        F1M2B=F1[i]
        F5M2B=F5[i]
        F6M2B=F6[i]
        F7M2B=F7[i]
        F8M2B=F8[i]
        F9M2B=F9[i]
        F10M2B=F10[i]
        F11M2B=F11[i]
        F12M2B=F12[i]
        F1M2B_li.append(F1M2B)
        F5M2B_li.append(F5M2B)
        F6M2B_li.append(F6M2B)
        F7M2B_li.append(F7M2B)
        F8M2B_li.append(F8M2B)
        F9M2B_li.append(F9M2B)
        F10M2B_li.append(F10M2B)
        F11M2B_li.append(F11M2B)
        F12M2B_li.append(F12M2B)
        
        
    if 10.0 < Mass[i] < 10.5 and Bar[i] == 0:
        F1M2NB=F1[i]
        F5M2NB=F5[i]
        F6M2NB=F6[i]
        F7M2NB=F7[i]
        F8M2NB=F8[i]
        F9M2NB=F9[i]
        F10M2NB=F10[i]
        F11M2NB=F11[i]
        F12M2NB=F12[i]
        F1M2NB_li.append(F1M2NB)
        F5M2NB_li.append(F5M2NB)
        F6M2NB_li.append(F6M2NB)
        F7M2NB_li.append(F7M2NB)
        F8M2NB_li.append(F8M2NB)
        F9M2NB_li.append(F9M2NB)
        F10M2NB_li.append(F10M2NB)
        F11M2NB_li.append(F11M2NB)
        F12M2NB_li.append(F12M2NB)
        
        
    if 10.5 < Mass[i] < 11.0 and Bar[i] == 1:
        F1M3B=F1[i]
        F5M3B=F5[i]
        F6M3B=F6[i]
        F7M3B=F7[i]
        F8M3B=F8[i]
        F9M3B=F9[i]
        F10M3B=F10[i]
        F11M3B=F11[i]
        F12M3B=F12[i]
        F1M3B_li.append(F1M3B)
        F5M3B_li.append(F5M3B)
        F6M3B_li.append(F6M3B)
        F7M3B_li.append(F7M3B)
        F8M3B_li.append(F8M3B)
        F9M3B_li.append(F9M3B)
        F10M3B_li.append(F10M3B)
        F11M3B_li.append(F11M3B)
        F12M3B_li.append(F12M3B)
        
        
    if 10.5 < Mass[i] < 11.0 and Bar[i] == 0:
        F1M3NB=F1[i]
        F5M3NB=F5[i]
        F6M3NB=F6[i]
        F7M3NB=F7[i]
        F8M3NB=F8[i]
        F9M3NB=F9[i]
        F10M3NB=F10[i]
        F11M3NB=F11[i]
        F12M3NB=F12[i]
        F1M3NB_li.append(F1M3NB)
        F5M3NB_li.append(F5M3NB)
        F6M3NB_li.append(F6M3NB)
        F7M3NB_li.append(F7M3NB)
        F8M3NB_li.append(F8M3NB)
        F9M3NB_li.append(F9M3NB)
        F10M3NB_li.append(F10M3NB)
        F11M3NB_li.append(F11M3NB)
        F12M3NB_li.append(F12M3NB)
        
    if 11.0 < Mass[i] < 11.5 and Bar[i] == 1:
        F1M4B=F1[i]
        F5M4B=F5[i]
        F6M4B=F6[i]
        F7M4B=F7[i]
        F8M4B=F8[i]
        F9M4B=F9[i]
        F10M4B=F10[i]
        F11M4B=F11[i]
        F12M4B=F12[i]
        F1M4B_li.append(F1M4B)
        F5M4B_li.append(F5M4B)
        F6M4B_li.append(F6M4B)
        F7M4B_li.append(F7M4B)
        F8M4B_li.append(F8M4B)
        F9M4B_li.append(F9M4B)
        F10M4B_li.append(F10M4B)
        F11M4B_li.append(F11M4B)
        F12M4B_li.append(F12M4B)
        
        
    if 11.0 < Mass[i] < 11.5 and Bar[i] == 0:
        F1M4NB=F1[i]
        F5M4NB=F5[i]
        F6M4NB=F6[i]
        F7M4NB=F7[i]
        F8M4NB=F8[i]
        F9M4NB=F9[i]
        F10M4NB=F10[i]
        F11M4NB=F11[i]
        F12M4NB=F12[i]
        F1M4NB_li.append(F1M4NB)
        F5M4NB_li.append(F5M4NB)
        F6M4NB_li.append(F6M4NB)
        F7M4NB_li.append(F7M4NB)
        F8M4NB_li.append(F8M4NB)
        F9M4NB_li.append(F9M4NB)
        F10M4NB_li.append(F10M4NB)
        F11M4NB_li.append(F11M4NB)
        F12M4NB_li.append(F12M4NB)
       
        
print('length M1_Barred:', len(F1M1B_li), len(F5M1B_li), len(F6M1B_li), len(F7M1B_li),len(F8M1B_li), len(F9M1B_li), len(F10M1B_li), len(F11M1B_li),len(F12M1B_li))
print('length M1_UNBarred:',len(F1M1NB_li), len(F5M1NB_li), len(F6M1NB_li), len(F7M1NB_li),len(F8M1NB_li), len(F9M1NB_li), len(F10M1NB_li), len(F11M1NB_li),len(F12M1NB_li))

print('length M2_Barred:', len(F1M2B_li), len(F5M2B_li), len(F6M2B_li), len(F7M2B_li),len(F8M2B_li), len(F9M2B_li), len(F10M2B_li), len(F11M2B_li),len(F12M2B_li))
print('length M2_UNBarred:',len(F1M2NB_li), len(F5M2NB_li), len(F6M2NB_li), len(F7M2NB_li),len(F8M2NB_li), len(F9M2NB_li), len(F10M2NB_li), len(F11M2NB_li),len(F12M2NB_li))

print('length M3_Barred:', len(F1M3B_li), len(F5M3B_li), len(F6M3B_li), len(F7M3B_li),len(F8M3B_li), len(F9M3B_li), len(F10M3B_li), len(F11M3B_li),len(F12M3B_li))
print('length M3_UNBarred:',len(F1M3NB_li), len(F5M3NB_li), len(F6M3NB_li), len(F7M3NB_li),len(F8M3NB_li), len(F9M3NB_li), len(F10M3NB_li), len(F11M3NB_li),len(F12M3NB_li))

print('length M4_Barred:', len(F1M4B_li), len(F5M4B_li), len(F6M4B_li), len(F7M4B_li),len(F8M4B_li), len(F9M4B_li), len(F10M4B_li), len(F11M4B_li),len(F12M4B_li))
print('length M4_UNBarred:',len(F1M4NB_li), len(F5M4NB_li), len(F6M4NB_li), len(F7M4NB_li),len(F8M4NB_li), len(F9M4NB_li), len(F10M4NB_li), len(F11M4NB_li),len(F12M4NB_li))


## y data for M1

y1M1B= statistics.mean(F1M1B_li)
y1M1NB= statistics.mean(F1M1NB_li)
y5M1B= statistics.mean(F5M1B_li)
y5M1NB= statistics.mean(F5M1NB_li)
y6M1B= statistics.mean(F6M1B_li)
y6M1NB= statistics.mean(F6M1NB_li)
y7M1B= statistics.mean(F7M1B_li)
y7M1NB= statistics.mean(F7M1NB_li)
y8M1B= statistics.mean(F8M1B_li)
y8M1NB= statistics.mean(F8M1NB_li)
y9M1B= statistics.mean(F9M1B_li)
y9M1NB= statistics.mean(F9M1NB_li)
y10M1B= statistics.mean(F10M1B_li)
y10M1NB= statistics.mean(F10M1NB_li)
y11M1B= statistics.mean(F11M1B_li)
y11M1NB= statistics.mean(F11M1NB_li)
y12M1B= statistics.mean(F12M1B_li)
y12M1NB= statistics.mean(F12M1NB_li)


print('y1M1b_nb',y1M1B, y1M1NB)
print('y5M1b_nb',y5M1B, y5M1NB)
print('y6M1b_nb',y6M1B, y6M1NB)
print('y7M1b_nb',y7M1B, y7M1NB)
print('y8M1b_nb',y8M1B, y8M1NB)
print('y9M1b_nb',y9M1B, y9M1NB)
print('y10M1b_nb',y10M1B, y10M1NB)
print('y11M1b_nb',y11M1B, y11M1NB)
print('y12M1b_nb',y12M1B, y12M1NB)


## y data for M2

y1M2B= statistics.mean(F1M2B_li)
y1M2NB= statistics.mean(F1M2NB_li)
y5M2B= statistics.mean(F5M2B_li)
y5M2NB= statistics.mean(F5M2NB_li)
y6M2B= statistics.mean(F6M2B_li)
y6M2NB= statistics.mean(F6M2NB_li)
y7M2B= statistics.mean(F7M2B_li)
y7M2NB= statistics.mean(F7M2NB_li)
y8M2B= statistics.mean(F8M2B_li)
y8M2NB= statistics.mean(F8M2NB_li)
y9M2B= statistics.mean(F9M2B_li)
y9M2NB= statistics.mean(F9M2NB_li)
y10M2B= statistics.mean(F10M2B_li)
y10M2NB= statistics.mean(F10M2NB_li)
y11M2B= statistics.mean(F11M2B_li)
y11M2NB= statistics.mean(F11M2NB_li)
y12M2B= statistics.mean(F12M2B_li)
y12M2NB= statistics.mean(F12M2NB_li)

print('y1M2b_nb',y1M2B, y1M2NB)
print('y5M2b_nb',y5M2B, y5M2NB)
print('y6M2b_nb',y6M2B, y6M2NB)
print('y7M2b_nb',y7M2B, y7M2NB)
print('y8M2b_nb',y8M2B, y8M2NB)
print('y9M2b_nb',y9M2B, y9M2NB)
print('y10M2b_nb',y10M2B, y10M2NB)
print('y11M2b_nb',y11M2B, y11M2NB)
print('y12M2b_nb',y12M2B, y12M2NB)

## y data for M3

y1M3B= statistics.mean(F1M3B_li)
y1M3NB= statistics.mean(F1M3NB_li)
y5M3B= statistics.mean(F5M3B_li)
y5M3NB= statistics.mean(F5M3NB_li)
y6M3B= statistics.mean(F6M3B_li)
y6M3NB= statistics.mean(F6M3NB_li)
y7M3B= statistics.mean(F7M3B_li)
y7M3NB= statistics.mean(F7M3NB_li)
y8M3B= statistics.mean(F8M3B_li)
y8M3NB= statistics.mean(F8M3NB_li)
y9M3B= statistics.mean(F9M3B_li)
y9M3NB= statistics.mean(F9M3NB_li)
y10M3B= statistics.mean(F10M3B_li)
y10M3NB= statistics.mean(F10M3NB_li)
y11M3B= statistics.mean(F11M3B_li)
y11M3NB= statistics.mean(F11M3NB_li)
y12M3B= statistics.mean(F12M3B_li)
y12M3NB= statistics.mean(F12M3NB_li)


print('y1M3b_nb',y1M3B, y1M3NB)
print('y5M3b_nb',y5M3B, y5M3NB)
print('y6M3b_nb',y6M3B, y6M3NB)
print('y7M3b_nb',y7M3B, y7M3NB)
print('y8M3b_nb',y8M3B, y8M3NB)
print('y9M3b_nb',y9M3B, y9M3NB)
print('y10M3b_nb',y10M3B, y10M3NB)
print('y11M3b_nb',y11M3B, y11M3NB)
print('y12M3b_nb',y12M3B, y12M3NB)

## y data for M4

y1M4B= statistics.mean(F1M4B_li)
y1M4NB= statistics.mean(F1M4NB_li)
y5M4B= statistics.mean(F5M4B_li)
y5M4NB= statistics.mean(F5M4NB_li)
y6M4B= statistics.mean(F6M4B_li)
y6M4NB= statistics.mean(F6M4NB_li)
y7M4B= statistics.mean(F7M4B_li)
y7M4NB= statistics.mean(F7M4NB_li)
y8M4B= statistics.mean(F8M4B_li)
y8M4NB= statistics.mean(F8M4NB_li)
y9M4B= statistics.mean(F9M4B_li)
y9M4NB= statistics.mean(F9M4NB_li)
y10M4B= statistics.mean(F10M4B_li)
y10M4NB= statistics.mean(F10M4NB_li)
y11M4B= statistics.mean(F11M4B_li)
y11M4NB= statistics.mean(F11M4NB_li)
y12M4B= statistics.mean(F12M4B_li)
y12M4NB= statistics.mean(F12M4NB_li)

print('y1M4b_nb',y1M4B, y1M4NB)
print('y5M4b_nb',y5M4B, y5M4NB)
print('y6M4b_nb',y6M4B, y6M4NB)
print('y7M4b_nb',y7M4B, y7M4NB)
print('y8M4b_nb',y8M4B, y8M4NB)
print('y9M4b_nb',y9M4B, y9M4NB)
print('y10M4b_nb',y10M4B, y10M4NB)
print('y11M4b_nb',y11M4B, y11M4NB)
print('y12M4b_nb',y12M4B, y12M4NB)


#sys.exit()

plt.figure(1)
fig, axs = plt.subplots(figsize=(8,3)) 

gs = gridspec.GridSpec(1,4)

ax1 = plt.subplot(gs[:,0]) 
x1= '7.0'
x2= '7.30'
x3='7.5'
x4='7.76'
x5='8.0'    
x6= '8.31'  
x7= '8.5'    
x8= '8.76'
x9='8.85'
x10='9.01'
x11='9.25'
x12='9.48'
x13='9.55'
x14='9.76'
x15='9.85'
x16='10.0'
x17='10.10'
x18='10.15' 
y1=100*(y1M1B) 
y2=100*(y1M1NB) 
y3=100*(y5M1B)
y4=100*(y5M1NB)
y5=100*(y6M1B)
y6=100*(y6M1NB)
y7=100*(y7M1B)
y8=100*(y7M1NB)
y9=100*(y8M1B) 
y10=100*(y8M1NB)
y11=100*(y9M1B)
y12=100*(y9M1NB)
y13=100*(y10M1B)
y14=100*(y10M1NB)
y15=100*(y11M1B)
y16=100*(y11M1NB)
y17=100*(y12M1B)
y18=100*(y12M1NB)



plt.bar(x1, y1, color='blue',width=0.5, label='Barred')
plt.bar(x2, y2, color='grey',width=0.5, label='Non-barred') 
plt.bar(x3, y3, color='blue',width=0.5)
plt.bar(x4,y4 ,color='grey',width=0.5)
plt.bar(x5,y5, color='blue', width=0.5)
plt.bar(x6,y6, color='grey', width=0.5) 
plt.bar(x7,y7, color='blue', width=0.5)
plt.bar(x8,y8, color='grey', width=0.5)
plt.bar(x9, y9, color='blue',width=0.5)
plt.bar(x10,y10 ,color='grey',width=0.5)
plt.bar(x11,y11, color='blue', width=0.5)
plt.bar(x12,y12, color='grey', width=0.5) 
plt.bar(x13,y13, color='blue', width=0.5)
plt.bar(x14,y14, color='grey', width=0.5)
plt.bar(x15,y15, color='blue', width=0.5)
plt.bar(x16,y16, color='grey', width=0.5) 
plt.bar(x17,y17, color='blue', width=0.5)
plt.bar(x18,y18, color='grey', width=0.5)
#plt.xlim(-3,10)
plt.ylim(0,60)
#plt.xticks(np.arange(-1, 7, step=1),fontsize=4)
plt.xticks(fontsize=4)
plt.xticks(rotation='vertical')
#y_vals = ax1.get_yticks()
#ax1.set_yticklabels(['{:1.0f}%'.format(log(x)) for x in y_vals])
plt.title("9.5 < M < 10.0",fontsize=6 )
plt.text(0.5,36,'2286(B)/',color='blue',fontsize=5)
plt.text(5.3,36,'7785(NB)',color='black',fontsize=5)
plt.xlabel('log(age/yr)') 
plt.ylabel('mass fraction[$\%$]')
plt.legend(fontsize=5, loc='upper left')

   
ax2 = plt.subplot(gs[:,1], sharey=ax1) 
x1= '7.0'
x2= '7.30'
x3='7.5'
x4='7.76'
x5='8.0'    
x6= '8.31'  
x7= '8.5'    
x8= '8.76'
x9='8.85'
x10='9.01'
x11='9.25'
x12='9.48'
x13='9.55'
x14='9.76'
x15='9.85'
x16='10.0'
x17='10.10'
x18='10.15' 
y1=100*(y1M2B) 
y2=100*(y1M2NB) 
y3=100*(y5M2B)
y4=100*(y5M2NB)
y5=100*(y6M2B)
y6=100*(y6M2NB)
y7=100*(y7M2B)
y8=100*(y7M2NB)
y9=100*(y8M2B) 
y10=100*(y8M2NB)
y11=100*(y9M2B)
y12=100*(y9M2NB)
y13=100*(y10M2B)
y14=100*(y10M2NB)
y15=100*(y11M2B)
y16=100*(y11M2NB)
y17=100*(y12M2B)
y18=100*(y12M2NB)


plt.bar(x1, y1, color='blue',width=0.5, label='Barred')
plt.bar(x2, y2, color='grey',width=0.5, label='Non-barred') 
plt.bar(x3, y3, color='blue',width=0.5)
plt.bar(x4,y4 ,color='grey',width=0.5)
plt.bar(x5,y5, color='blue', width=0.5)
plt.bar(x6,y6, color='grey', width=0.5) 
plt.bar(x7,y7, color='blue', width=0.5)
plt.bar(x8,y8, color='grey', width=0.5)
plt.bar(x9, y9, color='blue',width=0.5)
plt.bar(x10,y10 ,color='grey',width=0.5)
plt.bar(x11,y11, color='blue', width=0.5)
plt.bar(x12,y12, color='grey', width=0.5) 
plt.bar(x13,y13, color='blue', width=0.5)
plt.bar(x14,y14, color='grey', width=0.5)
plt.bar(x15,y15, color='blue', width=0.5)
plt.bar(x16,y16, color='grey', width=0.5) 
plt.bar(x17,y17, color='blue', width=0.5)
plt.bar(x18,y18, color='grey', width=0.5)
plt.xticks(fontsize=4,rotation='vertical')
#plt.ylim(0,1000)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.title("10.0 < M < 10.5",fontsize=6 )
plt.text(0.5,40,'2421(B)/',color='blue',fontsize=5)
plt.text(5.3,40,'6147(NB)',color='black',fontsize=5)
plt.xlabel('log(age/yr)', fontsize=10) 


ax3 = plt.subplot(gs[:,2], sharey=ax1) 
x1= '7.0'
x2= '7.30'
x3='7.5'
x4='7.76'
x5='8.0'    
x6= '8.31'  
x7= '8.5'    
x8= '8.76'
x9='8.85'
x10='9.01'
x11='9.25'
x12='9.48'
x13='9.55'
x14='9.76'
x15='9.85'
x16='10.0'
x17='10.10'
x18='10.15' 
y1=100*(y1M3B) 
y2=100*(y1M3NB) 
y3=100*(y5M3B)
y4=100*(y5M3NB)
y5=100*(y6M3B)
y6=100*(y6M3NB)
y7=100*(y7M3B)
y8=100*(y7M3NB)
y9=100*(y8M3B) 
y10=100*(y8M3NB)
y11=100*(y9M3B)
y12=100*(y9M3NB)
y13=100*(y10M3B)
y14=100*(y10M3NB)
y15=100*(y11M3B)
y16=100*(y11M3NB)
y17=100*(y12M3B)
y18=100*(y12M3NB)


plt.bar(x1, y1, color='blue',width=0.5, label='Barred')
plt.bar(x2, y2, color='grey',width=0.5, label='Non-barred') 
plt.bar(x3, y3, color='blue',width=0.5)
plt.bar(x4,y4 ,color='grey',width=0.5)
plt.bar(x5,y5, color='blue', width=0.5)
plt.bar(x6,y6, color='grey', width=0.5) 
plt.bar(x7,y7, color='blue', width=0.5)
plt.bar(x8,y8, color='grey', width=0.5)
plt.bar(x9, y9, color='blue',width=0.5)
plt.bar(x10,y10 ,color='grey',width=0.5)
plt.bar(x11,y11, color='blue', width=0.5)
plt.bar(x12,y12, color='grey', width=0.5) 
plt.bar(x13,y13, color='blue', width=0.5)
plt.bar(x14,y14, color='grey', width=0.5)
plt.bar(x15,y15, color='blue', width=0.5)
plt.bar(x16,y16, color='grey', width=0.5) 
plt.bar(x17,y17, color='blue', width=0.5)
plt.bar(x18,y18, color='grey', width=0.5)
plt.xticks( fontsize=4, rotation='vertical')
#plt.ylim(0,1000)
plt.setp(ax3.get_yticklabels(), visible=False)
plt.title("10.5 < M < 11.0",fontsize=6 )
plt.text(0.5,40,'1872(B)/',color='blue',fontsize=5)
plt.text(5.3,40,'3501(NB)',color='black',fontsize=5)
plt.xlabel('log(age/yr)')   

ax4 = plt.subplot(gs[:,3],sharey=ax1) 
x1= '7.0'
x2= '7.30'
x3='7.5'
x4='7.76'
x5='8.0'    
x6= '8.31'  
x7= '8.5'    
x8= '8.76'
x9='8.85'
x10='9.01'
x11='9.25'
x12='9.48'
x13='9.55'
x14='9.76'
x15='9.85'
x16='10.0'
x17='10.10'
x18='10.15' 
y1=100*(y1M4B) 
y2=100*(y1M4NB) 
y3=100*(y5M4B)
y4=100*(y5M4NB)
y5=100*(y6M4B)
y6=100*(y6M4NB)
y7=100*(y7M4B)
y8=100*(y7M4NB)
y9=100*(y8M4B) 
y10=100*(y8M4NB)
y11=100*(y9M4B)
y12=100*(y9M4NB)
y13=100*(y10M4B)
y14=100*(y10M4NB)
y15=100*(y11M4B)
y16=100*(y11M4NB)
y17=100*(y12M4B)
y18=100*(y12M4NB)

plt.bar(x1, y1, color='blue',width=0.5, label='Barred')
plt.bar(x2, y2, color='grey',width=0.5, label='Non-barred') 
plt.bar(x3, y3, color='blue',width=0.5)
plt.bar(x4,y4 ,color='grey',width=0.5)
plt.bar(x5,y5, color='blue', width=0.5)
plt.bar(x6,y6, color='grey', width=0.5) 
plt.bar(x7,y7, color='blue', width=0.5)
plt.bar(x8,y8, color='grey', width=0.5)
plt.bar(x9, y9, color='blue',width=0.5)
plt.bar(x10,y10 ,color='grey',width=0.5)
plt.bar(x11,y11, color='blue', width=0.5)
plt.bar(x12,y12, color='grey', width=0.5) 
plt.bar(x13,y13, color='blue', width=0.5)
plt.bar(x14,y14, color='grey', width=0.5)
plt.bar(x15,y15, color='blue', width=0.5)
plt.bar(x16,y16, color='grey', width=0.5) 
plt.bar(x17,y17, color='blue', width=0.5)
plt.bar(x18,y18, color='grey', width=0.5)
plt.xticks( fontsize=4, rotation='vertical')
#plt.ylim(0,1000)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.title("11.0 < M < 11.5",fontsize=6 )
plt.text(0.5,40,'432(B)/',color='blue',fontsize=5)
plt.text(5,40,'828(NB)',color='black',fontsize=5)
plt.xlabel('log(age/yr)')        


#plt.subplots_adjust(hspace=.0) ## remove the gap between the two plots
plt.subplots_adjust(left=0.2, bottom=0.18, top=0.85)

plt.savefig('sfr_12_massfraction.pdf')
plt.close()        
       
