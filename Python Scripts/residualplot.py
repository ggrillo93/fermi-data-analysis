#!/usr/bin/env python
from textwrap import wrap
import os
import pyfits
import numpy as np
import math
from pylab import *
from matplotlib import pyplot as plt
import csv
from subprocess import call
import glob

def fullText(text):
    with open(text) as f:
        lista=[]
        for line in f:
            lista.append(line)
    return lista

def columnToList(text,colnumber):
    with open(text) as f:
        lista=[]
        for line in f:
            spl=line.split()
            lista.append(spl[colnumber])
    return lista

def columnToList2(ptext,colnumber):
    lista=[]
    for line in ptext:
        spl=line.split()
        lista.append(spl[colnumber])
    return lista

def strToFloat(lista):
	nlista=[]
	for val in lista:
		try:
			newval=float(val)
		except ValueError:
			newval=0
		nlista.append(newval)
	return nlista
	
def positionLast (x,s):
    count = len(s)
    for i in s[::-1]:
        count -= 1
        if i == x:
            return count
    return None

def positionFirst(x,s):
    count=0
    for i in s:
        if i==x:
            return count
        else:
            count=count+1
    return None

def formatted(lista):
    dmlist=[]
    for dm in lista:
        newdm=''
        for ch in dm:
            if ch == '*':
                newdm=newdm+'e'
            elif ch == '^':
                newdm=newdm
            elif ch == '{':
                newdm=newdm
            elif ch ==',':
                newdm=newdm
            else:
                newdm=newdm+ch
        dmlist.append(newdm)
    return dmlist

def formatted2(lista):
    dmerrlist=[]
    for olderr in lista:
        new=olderr[:-1]
        dmerrlist.append(new)
    return dmerrlist

def findBeginning(tempofile):
    for n in range(0,len(tempofile)):
        if tempofile[n][:2]=='St':
            break
    return n

def epochDivider(epochloc,variable):
    varepochs=[]
    for pair in epochloc:
        startloc=pair[0]
        endloc=pair[1]
        varepoch=variable[startloc:endloc+1]
        varepochs.append(varepoch)
    return varepochs

def setMode(mode):
    mathinputdir="/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Input/"+pulsar+"/"+mode
    if not os.path.exists(mathinputdir):
        os.makedirs(mathinputdir)
    return mathinputdir

def findRMSres(tempofile):
	text=fullText(tempofile)
	for line in text:
		if line.startswith('RMS pre'):
			rmsres=line[-13:-1]
	if rmsres.startswith(' '):
		rmres=rmsres[1:]
	return rmsres+')'

def plotResiduals(par,tim,mode):
	# Run tempo2 on pulsar
	pulsar=par[:10]
	os.system('tempo2 -output general2 -f '+par+' '+tim+' -s "{freq} {post} {bat} {err}\n" > '+pulsar+'g2.txt')
	os.system('mv '+pulsar+'g2.txt /home/ggrillo93/Documents/Research/NANOGrav/DM/"iPython NB"/Input') # location of tempo output
	os.system('tempo2 -f '+par+' '+tim+' > '+pulsar+'reg.txt')
	os.system('mv '+pulsar+'reg.txt /home/ggrillo93/Documents/Research/NANOGrav/DM/"iPython NB"/Input')
	os.chdir("/home/ggrillo93/Documents/Research/NANOGrav/DM/iPython NB")
	
	# Open par file
	parfile="/home/ggrillo93/Documents/Research/NANOGrav/DM/Tempo2/Original/" + pulsar + "_NANOGrav_9yv1.gls.par"
	parameters=columnToList(parfile,0)
	parvalue=columnToList(parfile,1)
	
	if '-F' not in mode:
		# Locate lines containing lowest observation frequency for each epoch
		dmxf1loc=[]
		for n in range(0,len(parameters)):
			par=parameters[n]
			if par[:5] == 'DMXF1':
				dmxf1loc.append(n)
		
		# Locate lines containing highest observation frequency for each epoch
		dmxf2loc=[]
		for n in range(0,len(parameters)):
			par=parameters[n]
			if par[:5] == 'DMXF2':
				dmxf2loc.append(n)
		
		# Record lowest observation frequency for each epoch
		dmxf1=[]
		for loc in dmxf1loc:
			dmxfval=parvalue[loc]
			dmxf1.append(dmxfval)
		dmxf1=strToFloat(dmxf1)
		
		# Record highest observation frequency for each epoch
		dmxf2=[]
		for loc in dmxf2loc:
			dmxfval=parvalue[loc]
			dmxf2.append(dmxfval)
		dmxf2=strToFloat(dmxf2)    
		
		# Analyze dmxf1 and dmxf2 to generate list of observation frequencies
		f1min=min(dmxf1)
		obsfreqs=[]
		if f1min > 700 and f1min < 800:
			f1='820 MHz'
			obsfreqs.append(f1)
		if f1min < 500:
			f1='430 MHz'
			obsfreqs.append(f1)
		if f1min > 1000:
			f1='1410 MHz'
			obsfreqs.append(f1)
		f2max=max(dmxf2)
		if f2max > 1800 and f2max < 1900:
			f2='1500 MHz'
			obsfreqs.append(f2)
		if f2max > 1700 and f2max < 1800:
			f2='1410 MHz'
			obsfreqs.append(f2)
		if f2max > 2000:
			f2='2030 MHz'
			obsfreqs.append(f2)
		if f1 == '430 MHz' and f2=='2030 MHz':
			f2='1410 MHz'
			f3='2030 MHz'
			obsfreqs=[f1,f2,f3]
	
	# Open tempo2 output and delete unnecessary lines
	rawallbands="/home/ggrillo93/Documents/Research/NANOGrav/DM/iPython NB/Input/"+pulsar+"g2.txt"
	allbands=fullText(rawallbands)
	beg=findBeginning(allbands)
	allbands=allbands[beg+1:-2]
	
	# Extract frequencies, residuals, TOAs, and errors from tempo2 file
	freq=strToFloat(columnToList2(allbands,0))
	res=strToFloat(columnToList2(allbands,1))
	toas=strToFloat(columnToList2(allbands,2))
	err=strToFloat(columnToList2(allbands,3))
	data=[freq,res,toas,err]
	
	if '-F' not in mode:
		f1loc=[]
		f2loc=[]
		f3loc=[]
		for n in range(len(freq)):
			f=freq[n]
			if obsfreqs==['820 MHz','1500 MHz']:
				if f < 1000:
					f1loc.append(n)
				else:
					f2loc.append(n)
			if obsfreqs==['430 MHz','1410 MHz']:
				if f < 800:
					f1loc.append(n)
				else:
					f2loc.append(n)
			if obsfreqs==['1410 MHz','2030 MHz']:
				if f < 1761:
					f1loc.append(n)
				else:
					f2loc.append(n)
			if len(obsfreqs)==3:
				if f < 800:
					f1loc.append(n)
				elif f < 1761:
					f2loc.append(n)
				else:
					f3loc.append(n)
		if len(obsfreqs)==3:
			flocs=[f1loc,f2loc,f3loc]
		else:
			flocs=[f1loc,f2loc]
		
		# Arrange data based on frequency bands
		f1=[]
		f1res=[]
		f1toas=[]
		f1err=[]
		for loc in flocs[0]:
			f1.append(data[0][loc])
			f1res.append(data[1][loc])
			f1toas.append(data[2][loc])
			f1err.append(data[3][loc])
			
		f2=[]
		f2res=[]
		f2toas=[]
		f2err=[]
		for loc in flocs[1]:
			f2.append(data[0][loc])
			f2res.append(data[1][loc])
			f2toas.append(data[2][loc])
			f2err.append(data[3][loc])
			
		if len(obsfreqs)==3:
			f3=[]
			f3res=[]
			f3toas=[]
			f3err=[]
			for loc in flocs[2]:
				f3.append(data[0][loc])
				f3res.append(data[1][loc])
				f3toas.append(data[2][loc])
				f3err.append(data[3][loc])
			newdata=[[f1,f1res,f1toas,f1err],[f2,f2res,f2toas,f2err],[f3,f3res,f3toas,f3err]]
		else:
			newdata=[[f1,f1res,f1toas,f1err],[f2,f2res,f2toas,f2err]]
			
		regtempoout="/home/ggrillo93/Documents/Research/NANOGrav/DM/iPython NB/Input/"+pulsar+"reg.txt"
		colors=[]
		for freq in obsfreqs:
			if freq=='430 MHz':
				colors.append('green')
			if freq=='820 MHz':
				colors.append('red')
			if freq=='1410 MHz':
				colors.append('purple')
			if freq=='1500 MHz':
				colors.append('black')
			if freq=='2030 MHz':
				colors.append('orange')
		plt.figure()
		for r in range(len(newdata)):
			toas=newdata[r][2]
			res=np.multiply(newdata[r][1],10**6)
			err=newdata[r][3]
			plt.errorbar(toas,res, yerr=err, fmt='none',ecolor=colors[r])
			plt.scatter(toas,res,color=colors[r],marker='s',label=obsfreqs[r])
		plt.legend(fontsize=10,loc=1)
		plt.title(pulsar+ " Post-Fit "+mode+' (RMS residual '+findRMSres(regtempoout))
		plt.xlabel("Day (MJD)")
		plt.ylabel("Residuals (us)")
		plt.ylim(min(res)-10,max(res)+10)
		plt.minorticks_on
	
	else:
		regtempoout="/home/ggrillo93/Documents/Research/NANOGrav/DM/iPython NB/Input/"+pulsar+"reg.txt"
		plt.figure()
		toas=data[2]
		res=np.multiply(data[1],10**6)
		err=data[3]
		plt.errorbar(toas,res, yerr=err, fmt='none',ecolor='green')
		plt.scatter(toas,res,color='green',marker='s')
		if mode=='NNstar-F':
			plt.title(pulsar+ ' Post-Fit NN-F (RMS residual = '+findRMSres(regtempoout))	
		else:
			plt.title(pulsar+ " Post-Fit "+mode+' (RMS residual = '+findRMSres(regtempoout))
		plt.xlabel("Day (MJD)")
		plt.ylabel("Residuals (us)")
		plt.ylim(min(res)-10,max(res)+10)
		plt.minorticks_on
	plt.savefig('/home/ggrillo93/Documents/Research/NANOGrav/DM/Tempo2/Gamma/'+mode+'/Plots/'+pulsar+'.png')

# os.chdir('/home/ggrillo93/Documents/Research/NANOGrav/DM/Tempo2/Gamma/NN-NN/')
# pars=sorted(glob.glob('*mod.par'))
# tims=sorted(glob.glob('*9yv1.tim'))
# for n in range(len(pars)):
# 	os.chdir('/home/ggrillo93/Documents/Research/NANOGrav/DM/Tempo2/Gamma/NN-NN/')
# 	plotResiduals(pars[n],tims[n],'NN-NN')

os.chdir('/home/ggrillo93/Documents/Research/NANOGrav/DM/Tempo2/Gamma/NNstar-F/')
pars=sorted(glob.glob('*mod.par'))
tims=sorted(glob.glob('*chol.tim'))
for n in range(len(pars)):
	os.chdir('/home/ggrillo93/Documents/Research/NANOGrav/DM/Tempo2/Gamma/NNstar-F/')
	plotResiduals(pars[n],tims[n],'NNstar-F')
# 
# os.chdir('/home/ggrillo93/Documents/Research/NANOGrav/DM/Tempo2/Gamma/F-F/')
# pars=sorted(glob.glob('*mod.par'))
# tims=sorted(glob.glob('*chol.tim'))
# for n in range(len(pars)):
# 	os.chdir('/home/ggrillo93/Documents/Research/NANOGrav/DM/Tempo2/Gamma/F-F/')
# 	plotResiduals(pars[n],tims[n],'F-F')