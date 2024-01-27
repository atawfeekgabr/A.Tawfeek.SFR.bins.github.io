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





c1=20*10**6  #age bin at SFR1 =7.30
c2=570*10**6  ##=8.75
c3=7.5*10**9 ##=9.87
c4=13.5*10**9 ##10.13



S1M1B_li=[]
S1M1NB_li=[]

S1M2B_li=[]
S1M2NB_li=[]

S1M3B_li=[]
S1M3NB_li=[]

S1M4B_li=[]
S1M4NB_li=[]


S2M1B_li=[]
S3M1B_li=[]
S4M1B_li=[]

S2M2B_li=[]
S3M2B_li=[]
S4M2B_li=[]

S2M3B_li=[]
S3M3B_li=[]
S4M3B_li=[]

S2M4B_li=[]
S3M4B_li=[]
S4M4B_li=[]


S2M1NB_li=[]
S3M1NB_li=[]
S4M1NB_li=[]

S2M2NB_li=[]
S3M2NB_li=[]
S4M2NB_li=[]

S2M3NB_li=[]
S3M3NB_li=[]
S4M3NB_li=[]

S2M4NB_li=[]
S3M4NB_li=[]
S4M4NB_li=[]

F1M1B_li=[]
F1M1NB_li=[]
F1M2B_li=[]
F1M2NB_li=[]
F1M3B_li=[]
F1M3NB_li=[]
F1M4B_li=[]
F1M4NB_li=[]


F2M1B_li=[]
F2M1NB_li=[]
F2M2B_li=[]
F2M2NB_li=[]
F2M3B_li=[]
F2M3NB_li=[]
F2M4B_li=[]
F2M4NB_li=[]

F3M1B_li=[]
F3M1NB_li=[]
F3M2B_li=[]
F3M2NB_li=[]
F3M3B_li=[]
F3M3NB_li=[]
F3M4B_li=[]
F3M4NB_li=[]


F4M1B_li=[]
F4M1NB_li=[]
F4M2B_li=[]
F4M2NB_li=[]
F4M3B_li=[]
F4M3NB_li=[]
F4M4B_li=[]
F4M4NB_li=[]


d='disk_SFR_history_Mtype.dat'
data=np.loadtxt(d,dtype={'names':('WINGS', 'Weights','Bar','SFR1', 'SFR2','SFR3', 'SFR4','logM','M1','M2','M3','M4','Mtotal','F1','F2','F3','F4','Mtype'),'formats': ('O','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f')} ,skiprows=1)

Bar=data['Bar']
SFR1=data['SFR1']
SFR2=data['SFR2']
SFR3=data['SFR3']
SFR4=data['SFR4']
Mass=data['logM']
M1=data['M1']
M2=data['M2']
M3=data['M3']
M4=data['M4']
Mtotal=data['Mtotal']
F1=data['F1']
F2=data['F2']
F3=data['F3']
F4=data['F4']
Mtype=data['Mtype']




for i in range(len(Mtype)):
    if -5 < Mtype[i] < -2.5 and Bar[i] == 1:
        S1M1B= SFR1[i]
        S2M1B= SFR2[i]  
        S3M1B=SFR3[i]  
        S4M1B=SFR4[i] 
        S1M1NB= SFR1[i]
        S1M1NB_li.append(S1M1NB)
        S1M1B_li.append(S1M1B)
        S2M1B_li.append(S2M1B)
        S3M1B_li.append(S3M1B)
        S4M1B_li.append(S4M1B)
        F1M1B=F1[i]
        F2M1B=F2[i]
        F3M1B=F3[i]
        F4M1B=F4[i]
        F1M1B_li.append(F1M1B)
        F2M1B_li.append(F2M1B)
        F3M1B_li.append(F3M1B)
        F4M1B_li.append(F4M1B)

        
        
        
        
    if -5 < Mtype[i] < -2.5 and Bar[i] == 0:
        S1M1NB= SFR1[i]
        S2M1NB= SFR2[i] 
        S3M1NB=SFR3[i] 
        S4M1NB=SFR4[i] 
        S1M1NB_li.append(S1M1NB)
        S2M1NB_li.append(S2M1NB)
        S3M1NB_li.append(S3M1NB)
        S4M1NB_li.append(S4M1NB)
        F1M1NB=F1[i]
        F2M1NB=F2[i]
        F3M1NB=F3[i]
        F4M1NB=F4[i]
        F1M1NB_li.append(F1M1NB)
        F2M1NB_li.append(F2M1NB)
        F3M1NB_li.append(F3M1NB)
        F4M1NB_li.append(F4M1NB)
        
        
    if -2.5 < Mtype[i] < 0 and Bar[i] == 1:
        S2M2B= SFR2[i]
        S3M2B=SFR3[i]
        S4M2B=SFR4[i]
        S1M2B= SFR1[i]
        S1M2B_li.append(S1M2B)
        S2M2B_li.append(S2M2B)
        S3M2B_li.append(S3M2B)
        S4M2B_li.append(S4M2B)
        F1M2B=F1[i]
        F2M2B=F2[i]
        F3M2B=F3[i]
        F4M2B=F4[i]
        F1M2B_li.append(F1M2B)
        F2M2B_li.append(F2M2B)
        F3M2B_li.append(F3M2B)
        F4M2B_li.append(F4M2B)
        
        
    if -2.5 < Mtype[i] < 0 and Bar[i] == 0:
        S2M2NB= SFR2[i]
        S3M2NB=SFR3[i]
        S4M2NB=SFR4[i]
        S1M2NB= SFR1[i]
        S1M2NB_li.append(S1M2NB)
        S2M2NB_li.append(S2M2NB)
        S3M2NB_li.append(S3M2NB)
        S4M2NB_li.append(S4M2NB)
        F1M2NB=F1[i]
        F2M2NB=F2[i]
        F3M2NB=F3[i]
        F4M2NB=F4[i]
        F1M2NB_li.append(F1M2NB)
        F2M2NB_li.append(F2M2NB)
        F3M2NB_li.append(F3M2NB)
        F4M2NB_li.append(F4M2NB)
        
        
    if 0 < Mtype[i] < 2.5 and Bar[i] == 1:
        S2M3B= SFR2[i]
        S3M3B=SFR3[i]
        S4M3B=SFR4[i]
        S1M3B= SFR1[i]
        S1M3B_li.append(S1M3B)
        S2M3B_li.append(S2M3B)
        S3M3B_li.append(S3M3B)
        S4M3B_li.append(S4M3B)
        F1M3B=F1[i]
        F2M3B=F2[i]
        F3M3B=F3[i]
        F4M3B=F4[i]
        F1M3B_li.append(F1M3B)
        F2M3B_li.append(F2M3B)
        F3M3B_li.append(F3M3B)
        F4M3B_li.append(F4M3B)
        
        
    if 0 < Mtype[i] < 2.5 and Bar[i] == 0:
        S2M3NB= SFR2[i]
        S3M3NB=SFR3[i]
        S4M3NB=SFR4[i]
        S1M3NB= SFR1[i]
        S1M3NB_li.append(S1M3NB)
        S2M3NB_li.append(S2M3NB)
        S3M3NB_li.append(S3M3NB)
        S4M3NB_li.append(S4M3NB)
        F1M3NB=F1[i]
        F2M3NB=F2[i]
        F3M3NB=F3[i]
        F4M3NB=F4[i]
        F1M3NB_li.append(F1M3NB)
        F2M3NB_li.append(F2M3NB)
        F3M3NB_li.append(F3M3NB)
        F4M3NB_li.append(F4M3NB)
        
        
    if 2.5 < Mtype[i] < 6 and Bar[i] == 1:
        S2M4B= SFR2[i]
        S3M4B=SFR3[i]
        S4M4B=SFR4[i]
        S1M4B= SFR1[i]
        S1M4B_li.append(S1M4B)
        S2M4B_li.append(S2M4B)
        S3M4B_li.append(S3M4B)
        S4M4B_li.append(S4M4B)
        F1M4B=F1[i]
        F2M4B=F2[i]
        F3M4B=F3[i]
        F4M4B=F4[i]
        F1M4B_li.append(F1M4B)
        F2M4B_li.append(F2M4B)
        F3M4B_li.append(F3M4B)
        F4M4B_li.append(F4M4B)
        
        
    if 2.5 < Mtype[i] < 6 and Bar[i] == 0:
        S2M4NB= SFR2[i]
        S3M4NB=SFR3[i]
        S4M4NB=SFR4[i]
        S1M4NB= SFR1[i]
        S1M4NB_li.append(S1M4NB)
        S2M4NB_li.append(S2M4NB)
        S3M4NB_li.append(S3M4NB)
        S4M4NB_li.append(S4M4NB)
        F1M4NB=F1[i]
        F2M4NB=F2[i]
        F3M4NB=F3[i]
        F4M4NB=F4[i]
        F1M4NB_li.append(F1M4NB)
        F2M4NB_li.append(F2M4NB)
        F3M4NB_li.append(F3M4NB)
        F4M4NB_li.append(F4M4NB)
       
        
print('length M1_Barred:', len(S1M1B_li), len(S2M1B_li), len(S3M1B_li), len(S4M1B_li))
print('length M1_UNBarred:',len(S1M1NB_li),len(S2M1NB_li), len(S3M1NB_li), len(S4M1NB_li))

print('length M2_Barred:', len(S1M2B_li), len(S2M2B_li), len(S3M2B_li), len(S4M2B_li))
print('length M2_UNBarred:',len(S1M2NB_li),len(S2M2NB_li), len(S3M2NB_li), len(S4M2NB_li))

print('length M3_Barred:', len(S1M3B_li), len(S2M3B_li), len(S3M3B_li), len(S4M3B_li))
print('length M3_UNBarred:',len(S1M3NB_li),len(S2M3NB_li), len(S3M3NB_li), len(S4M3NB_li))

print('length M4_Barred:', len(S1M4B_li), len(S2M4B_li), len(S3M4B_li), len(S4M4B_li))
print('length M4_UNBarred:',len(S1M4NB_li),len(S2M4NB_li), len(S3M4NB_li), len(S4M4NB_li))


## y data for M1

y1M1B= statistics.mean(S1M1B_li)
y1M1NB= statistics.mean(S1M1NB_li)
y2M1B= statistics.mean(S2M1B_li)
y2M1NB= statistics.mean(S2M1NB_li)
y3M1B= statistics.mean(S3M1B_li)
y3M1NB= statistics.mean(S3M1NB_li)
y4M1B= statistics.mean(S4M1B_li)
y4M1NB= statistics.mean(S4M1NB_li)

print('y1M1b_nb',y1M1B, y1M1NB)
print('y2M1b_nb',y2M1B, y2M1NB)
print('y3M1b_nb',y3M1B, y3M1NB)
print('y4M1b_nb',y4M1B, y4M1NB)

## y data for M2

y1M2B= statistics.mean(S1M2B_li)
y1M2NB= statistics.mean(S1M2NB_li)
y2M2B= statistics.mean(S2M2B_li)
y2M2NB= statistics.mean(S2M2NB_li)
y3M2B= statistics.mean(S3M2B_li)
y3M2NB= statistics.mean(S3M2NB_li)
y4M2B= statistics.mean(S4M2B_li)
y4M2NB= statistics.mean(S4M2NB_li)

print('y1M2b_nb',y1M2B, y1M2NB)
print('y2M2b_nb',y2M2B, y2M2NB)
print('y3M2b_nb',y3M2B, y3M2NB)
print('y4M2b_nb',y4M2B, y4M2NB)



## y data for M3

y1M3B= statistics.mean(S1M3B_li)
y1M3NB= statistics.mean(S1M3NB_li)
y2M3B= statistics.mean(S2M3B_li)
y2M3NB= statistics.mean(S2M3NB_li)
y3M3B= statistics.mean(S3M3B_li)
y3M3NB= statistics.mean(S3M3NB_li)
y4M3B= statistics.mean(S4M3B_li)
y4M3NB= statistics.mean(S4M3NB_li)

print('y1M3b_nb',y1M3B, y1M3NB)
print('y2M3b_nb',y2M3B, y2M3NB)
print('y3M3b_nb',y3M3B, y3M3NB)
print('y4M3b_nb',y4M3B, y4M3NB)

## y data for M4

y1M4B= statistics.mean(S1M4B_li)
y1M4NB= statistics.mean(S1M4NB_li)
y2M4B= statistics.mean(S2M4B_li)
y2M4NB= statistics.mean(S2M4NB_li)
y3M4B= statistics.mean(S3M4B_li)
y3M4NB= statistics.mean(S3M4NB_li)
y4M4B= statistics.mean(S4M4B_li)
y4M4NB= statistics.mean(S4M4NB_li)

print('y1M4b_nb',y1M4B, y1M4NB)
print('y2M4b_nb',y2M4B, y2M4NB)
print('y3M4b_nb',y3M4B, y3M4NB)
print('y4M4b_nb',y4M4B, y4M4NB)


#sys.exit()

plt.figure(1)
fig, axs = plt.subplots(figsize=(6,3)) 

gs = gridspec.GridSpec(1,4)

ax1 = plt.subplot(gs[:,0]) 
x1= '7.0'
x2= '7.30'
x3='8.50'
x4='8.75'
x5='9.6'    #(S3M1B_li)
x6= '9.87'  #(S3M1NB_li)
x7= '10'    # (S4M1B_li)
x8= '10.13' # (S4M1NB_li)
y1=y1M1B 
y2=y1M1NB 
y3=y2M1B
y4=y2M1NB
y5=y3M1B
y6=y3M1NB
y7=y4M1B
y8=y4M1NB
plt.bar(x1, y1, color='blue',width=0.5, label='Barred')
plt.bar(x2, y2, color='grey',width=0.5, label='non-barred') 
plt.bar(x3, y3, color='blue',width=0.5)
plt.bar(x4,y4 ,color='grey',width=0.5)
plt.bar(x5,y5, color='blue', width=0.5)
plt.bar(x6,y6, color='grey', width=0.5) 
plt.bar(x7,y7, color='blue', width=0.5)
plt.bar(x8,y8, color='grey', width=0.5)
#plt.xlim(-3,10)
plt.ylim(0,2.0)
#plt.xticks(np.arange(-1, 7, step=1),fontsize=4)
plt.xticks(fontsize=4)
plt.xticks(rotation='vertical')
#y_vals = ax1.get_yticks()
#ax1.set_yticklabels(['{:1.0f}%'.format(x * 100) for x in y_vals])
plt.title("-5 < Mt < -2.5",fontsize=6 )
#plt.text(4,2.9,'922(B)/',color='blue',fontsize=5)
#plt.text(4,2.7,'2975(NB)',color='black',fontsize=5)
plt.xlabel('log(age/yr)') 
plt.ylabel('SFR')
plt.legend(fontsize=5, loc='upper left')

   
ax2 = plt.subplot(gs[:,1], sharey=ax1) 
x1= '7.0'
x2= '7.30'
x3='8.50'
x4='8.75'
x5='9.6'    #(S3M2B_li)
x6= '9.87'  #(S3M2NB_li)
x7= '10'    # (S4M2B_li)
x8= '10.13' # (S4M2NB_li)
y1=y1M2B 
y2=y1M2NB 
y3=y2M2B 
y4=y2M2NB 
y5=y3M2B 
y6=y3M2NB 
y7=y4M2B 
y8=y4M2NB 
plt.bar(x1, y1, color='blue',width=0.5, label='S1B')
plt.bar(x2, y2, color='grey',width=0.5, label='S1NB') 
plt.bar(x3, y3, color='blue',width=0.5, label='S2B')
plt.bar(x4,y4, color='grey',width=0.5, label='S2NB')
plt.bar(x5,y5, color='blue', width=0.5, label='S3B')
plt.bar(x6,y6, color='grey', width=0.5, label='S3NB') 
plt.bar(x7,y7, color='blue', width=0.5, label='S4B')
plt.bar(x8,y8, color='grey', width=0.5,label='S4NB')
plt.xticks(fontsize=4,rotation='vertical')
#plt.ylim(0,1000)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.title("-2.5 < Mt < 0",fontsize=6 )
#plt.text(0.5,3,'930(B)/',color='blue',fontsize=5)
#plt.text(3,3,'2295(NB)',color='black',fontsize=5)
plt.xlabel('log(age/yr)', fontsize=10) 


ax3 = plt.subplot(gs[:,2], sharey=ax1) 
x1= '7.0'
x2= '7.30'
x3='8.50'
x4='8.75'
x5='9.6'    #(S3M3B_li)
x6= '9.87'  #(S3M3NB_li)
x7= '10'    # (S4M3B_li)
x8= '10.13' # (S4M3NB_li)
y1=y1M3B 
y2=y1M3NB 
y3=y2M3B 
y4=y2M3NB 
y5=y3M3B 
y6=y3M3NB 
y7=y4M3B 
y8=y4M3NB 
plt.bar(x1, y1, color='blue',width=0.5, label='S1B')
plt.bar(x2, y2, color='grey',width=0.5, label='S1NB') 
plt.bar(x3, y3, color='blue',width=0.5, label='S2B')
plt.bar(x4,y4, color='grey',width=0.5, label='S2NB')
plt.bar(x5,y5, color='blue', width=0.5, label='S3B')
plt.bar(x6,y6, color='grey', width=0.5, label='S3NB') 
plt.bar(x7,y7, color='blue', width=0.5, label='S4B')
plt.bar(x8,y8, color='grey', width=0.5,label='S4NB')
plt.xticks( fontsize=4, rotation='vertical')
#plt.ylim(0,1000)
plt.setp(ax3.get_yticklabels(), visible=False)
plt.title("0 < Mt < 2.5",fontsize=6 )
#plt.text(0.5,3,'690(B)/',color='blue',fontsize=5)
#plt.text(3,3,'1258(NB)',color='black',fontsize=5)
plt.xlabel('log(age/yr)')   

ax4 = plt.subplot(gs[:,3],sharey=ax1) 
x1= '7.0'
x2= '7.30'
x3='8.50'
x4='8.75'
x5='9.6'    #(S3M4B_li)
x6= '9.87'  #(S3M4NB_li)
x7= '10'    # (S4M4B_li)
x8= '10.13' # (S4M4NB_li)
y1=y1M4B 
y2=y1M4NB 
y3=y2M4B 
y4=y2M4NB 
y5=y3M4B 
y6=y3M4NB 
y7=y4M4B 
y8=y4M4NB 
plt.bar(x1, y1, color='blue',width=0.5, label='S1B')
plt.bar(x2, y2, color='grey',width=0.5, label='S1NB') 
plt.bar(x3, y3, color='blue',width=0.5, label='S2B')
plt.bar(x4,y4, color='grey',width=0.5, label='S2NB')
plt.bar(x5,y5, color='blue', width=0.5, label='S3B')
plt.bar(x6,y6, color='grey', width=0.5, label='S3NB') 
plt.bar(x7,y7, color='blue', width=0.5, label='S4B')
plt.bar(x8,y8, color='grey', width=0.5,label='S4NB')
plt.xticks( fontsize=4, rotation='vertical')
#plt.ylim(0,1000)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.title("2.5 < Mt < 6",fontsize=6 )
#plt.text(0.5,3,'158(B)/',color='blue',fontsize=5)
#plt.text(3,3,'295(NB)',color='black',fontsize=5)
plt.xlabel('log(age/yr)')        


#plt.subplots_adjust(hspace=.0) ## remove the gap between the two plots
plt.subplots_adjust(left=0.2, bottom=0.18, top=0.85)

plt.savefig('SFR_4_Mtype.pdf')
plt.close()        
       
