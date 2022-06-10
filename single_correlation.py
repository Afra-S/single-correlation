import numpy as np
import sys
import subprocess
from math import *
from scipy.stats import circmean
from scipy.stats import circstd

from scipy.stats import pearsonr
from scipy.stats import kendalltau
from scipy.stats import spearmanr
from astropy.stats import circcorrcoef
from astropy import units as u

from scipy.spatial import distance

logdata1=input("First set: MD data or experimental data? [MD/exp] ")
logang1=False
logang2=False
if(logdata1.lower()=='md'):
    name1=input("Enter name file MD ")
    tet=input("Are you considering angular averages?[Y/N] ")
    if(tet.lower()=='y'):
        logang1=True
else:
    name1=input("Enter name file exp data ")

logdata2=input("Second set: MD data or experimental data? [MD/exp] ")
if(logdata2.lower()=='md'):
    name2=input("Enter name file MD ")
    tet=input("Are you considering angular averages?[Y/N] ")
    if(tet.lower()=='y'):
        logang2=True
else:
    name2=input("Enter name file exp data ")
logmask='n'
##Correlation with a mask
if((logdata1=='exp' and logdata2=='md') or (logdata2=='exp' and logdata1=='md')):
    logmask=input("Do you want to use discrete levels for the correlation?[Y/N] ")
    if(logmask.lower()=='n'):
        if(logang1==False and logang2==False):
            corrtype=input("What type of correlations do you want to compute ( Pearsons [P], Kendall [K] or Spearman [S])? [P/K/S] ")
    else:
        nlevel=int(input("How many levels do you want to use (maximum 3 level)? "))
        leveln=np.zeros((nlevel))

        for ll in range(0,nlevel):
            leveln[ll]=float(input("Value "+str(ll+1)+'-th level'))
            
corrtype=input("What type of correlations do you want to compute ( Pearsons [P], Kendall [K] or Spearman [S])? [P/K/S] ")
    
file1=open(name1,'r')
h=1000
data1=np.zeros((h))
ki=-1
for line in file1:
    ki=ki+1
    #print(line)
    data1[ki]=float(line)
file2=open(name2,'r')
h=1000
data2=np.zeros((h))
ki=-1
for line in file2:
    ki=ki+1
    #print(line)
    data2[ki]=float(line)
    
logcorrtype=False
ktot=ki+1
xdata=np.zeros((ktot))
ydata=np.zeros((ktot))
if(logmask.lower()=='n'):
    if((logdata1.lower()=='md') and (logdata2.lower()=='md')):
        if(logang1==True and logang2==True):
            corrvalue=circcorrcoef(data1[0:ki+1]*u.deg, data2[0:ki+1]*u.deg)
        elif((logang1==True and logang2==False)):
            ydata=np.arctan(data2[0:ki+1])
            corrvalue=circcorrcoef(data1[0:ki+1]*u.deg,ydata)
        elif ((logang1==False and logang2==True)):
            xdata=np.arctan(data1[0:ki+1])
            corrvalue=circcorrcoef(xdata, data2[0:ki+1]*u.deg)
        else:
            logcorrtype=True
            xdata1=data1[0:ki+1]
            ydata1=data2[0:ki+1]
    elif((logdata1.lower()=='md') and (logdata2.lower()=='exp')):
        xdata=np.zeros((ktot))
        ydata=np.zeros((ktot))
        ti=-1
        for ll in range(0,ktot):
            if(data2[ll] >-0.00001):
                ti=ti+1
                xdata[ti]=data1[ll]
                ydata[ti]=data2[ll]
        if(logang1==True):
            xdata1=xdata[0:ti+1]*u.deg
            ydata1=np.arctan(ydata[0:ti+1])
            corrvalue=circcorrcoef(xdata1,ydata1)
        else:
            logcorrtype=True
            xdata1=xdata[0:ti+1]
            ydata1=ydata[0:ti+1]
    elif((logdata2.lower()=='md') and (logdata1.lower()=='exp')):
        xdata=np.zeros((ktot))
        ydata=np.zeros((ktot))
        ti=-1
        for ll in range(0,ktot):
            if(data1[ll] >-0.00001):
                ti=ti+1
                xdata[ti]=data1[ll]
                ydata[ti]=data2[ll]
        if(logang2==True):
            ydata1=ydata[0:ti+1]*u.deg
            xdata1=np.arctan(xdata[0:ti+1])
            corrvalue=circcorrcoef(xdata1,ydata1)
        else:
            logcorrtype=True
            xdata1=xdata[0:ti+1]
            ydata1=ydata[0:ti+1]
    else:
        ti=-1
        for ll in range(0,ktot):
            if(data1[ll] >-0.00001 and data2[ll] >-0.00001):
                ti=ti+1
                xdata[ti]=data1[ll]
                ydata[ti]=data2[ll]
        logcorrtype=True
        xdata1=xdata[0:ti+1]
        ydata1=ydata[0:ti+1]


    if(logcorrtype==True):
        if(corrtype.lower()=='p'):
            corrvalue,p_value=pearsonr(xdata1,ydata1)
            print('Pearson correlation value= ',corrvalue)
            print('p-value= ',p_value)
        elif(corrtype.lower()=='k'):
            corrvalue,p_value=kendalltau(xdata1,ydata1)
            print('Kendall tau coefficient = ',corrvalue)
            print('p-value= ',p_value)
        elif(corrtype.lower()=='s'):
            corrvalue,p_value=spearmanr(xdata1,ydata1)
            print('Spearman correlation value= ',corrvalue)
            print('p-value= ',p_value)
    else:
        print(corrvalue)
    
else:
    levelnexp=np.zeros((nlevel))
    if(nlevel==2):
        levelnexp[0]=0
        levelnexp[1]=0.4
    else:
        levelnexp[2]=0.7
    if(logdata1.lower()=='exp'):
        xdata=np.zeros((ktot))
        ydata=np.zeros((ktot))
        ti=-1
        
        for ll in range(0,ktot):
            if(data1[ll] >-0.00001):
                ti=ti+1
                logd1=False
                for ttt in range(1,nlevel):
                    if(data1[ll] < levelnexp[ttt]):
                        xdata[ti]=ttt-1
                        logd1=True
                        break
                if(logd1==False):
                    xdata[ti]=nlevel-1
                logd1=False
                for ttt in range(1,nlevel):
                    if(data2[ll] < leveln[ttt]):
                        ydata[ti]=ttt-1
                        logd1=True
                        break
                if(logd1==False):
                    ydata[ti]=nlevel-1
                print(ti,data2[ll],xdata[ti],ydata[ti])
    xdata1=xdata[0:ti+1]
    ydata1=ydata[0:ti+1]
    if(corrtype.lower()=='p'):
        corrvalue,p_value=pearsonr(xdata1,ydata1)
        print('Pearson correlation value= ',corrvalue)
        print('p-value= ',p_value)
    elif(corrtype.lower()=='k'):
        corrvalue,p_value=kendalltau(xdata1,ydata1)
        print('Kendall tau coefficient = ',corrvalue)
        print('p-value= ',p_value)
    elif(corrtype.lower()=='s'):
        corrvalue,p_value=spearmanr(xdata1,ydata1)
        print('Spearman correlation value= ',corrvalue)
        print('p-value= ',p_value)

    dist1=distance.euclidean(xdata1,ydata1)
    print('Euclidean distance = ',dist1)
    dist2=sum(abs(aaa-bbb) for aaa,bbb in zip(xdata1,ydata1))
    print('Manhattan distance = ',dist2)
    
