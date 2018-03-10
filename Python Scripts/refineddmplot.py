#!/usr/bin/env python
import os
import pyfits
import numpy as np
import math
from pylab import *
from matplotlib import pyplot as plt
import csv
from subprocess import call

fullpulsar=raw_input('Enter pulsar name: ')
modes=raw_input('Enter modes (Crossband, Dualband and/or Inband): ')
pulsar=fullpulsar[:5]

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

# Run tempo2 on pulsar
os.chdir("/home/ggrillo93/Documents/Research/NANOGrav/DM/Tempo2/Original/") # tempo files location
first='tempo2 -output general2 -f '+fullpulsar+'_NANOGrav_9yv1.gls.par '+fullpulsar+'_NANOGrav_9yv1.tim '
second='-s "{freq} {pre} {bat} {err}\n" > '+pulsar+'.txt'
os.system(first+second)
os.system('mv '+pulsar+'.txt /home/ggrillo93/Documents/Research/NANOGrav/DM/"iPython NB"/Input') # location of tempo output
os.chdir("/home/ggrillo93/Documents/Research/NANOGrav/DM/iPython NB")

# Open par file
parfile="/home/ggrillo93/Documents/Research/NANOGrav/DM/Tempo2/Original/" + fullpulsar + "_NANOGrav_9yv1.gls.par"
parameters=columnToList(parfile,0)
parvalue=columnToList(parfile,1)

dmxfile=fullText('/home/ggrillo93/Documents/Research/NANOGrav/DM/Tempo2/Gamma/NNstar-F2/'+fullpulsar+'/'+fullpulsar+'.dmt')[3:]

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

# Locate epochs containing observations in one band only
diff=np.subtract(dmxf2,dmxf1)
singlebandloc=[]
for n in range(0,len(diff)):
    if diff[n] < 1000:
        singlebandloc.append(n)

# Locate lines containing DMX values for each epoch
dmxloc=[]
for n in range(0,len(parameters)):
    par=parameters[n]
    if par[:5] == 'DMX_0':
        dmxloc.append(n)

# Record DMX value in tempo2 format
tdmx=[]
for loc in dmxloc:
    dmxval=parvalue[loc]
    tdmx.append(dmxval)

# Transform DMX values from par file format to Python
dmx=[]
for dmxvalue in tdmx:
    newdmxvalue=''
    for ch in dmxvalue:
        if ch != 'D':
            newdmxvalue=newdmxvalue+ch
        else:
            newdmxvalue=newdmxvalue+'e'
    dmx.append(newdmxvalue)
dmx=strToFloat(dmx)

# Locate lines containing epoch ranges
dmxrangeloc=[]
for n in range(0,len(parameters)):
    par=parameters[n]
    if par[:4] == 'DMXR':
        dmxrangeloc.append(n)

# Record ranges values for each epoch
dmxr=[]
for loc in dmxrangeloc:
    day=parvalue[loc]
    dmxr.append(day)
dmxr=strToFloat(dmxr)

# Rearrange ranges values in pairs
dmxranges=[]
for n in range(0,len(dmxr)-1,2):
    r=[dmxr[n],dmxr[n+1]]
    dmxranges.append(r)

# Assign a particular MJD to each epoch
epochdays=[]
for dmxpair in dmxranges:
    av=(dmxpair[0]+dmxpair[1])/2.0
    epochdays.append(av)

# Modify lists to remove data for epochs containing only single band observations
adjsbloc=[]
for n in range(0,len(singlebandloc)):
    new=singlebandloc[n]-n
    adjsbloc.append(new)
    
for loc in adjsbloc:
	del epochdays[loc]
	del dmxranges[loc]
	del dmxf1[loc]
	del dmxf2[loc]
	del dmx[loc]
	del dmxfile[loc]

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
rawallbands="/home/ggrillo93/Documents/Research/NANOGrav/DM/iPython NB/Input/"+pulsar+".txt"
allbands=fullText(rawallbands)
beg=findBeginning(allbands)
allbands=allbands[beg+1:-2]

# Extract frequencies, residuals, TOAs, and errors from tempo2 file
freq=strToFloat(columnToList2(allbands,0))
res=strToFloat(columnToList2(allbands,1))
toas=strToFloat(columnToList2(allbands,2))
err=strToFloat(columnToList2(allbands,3))

# Sort quantities based on TOAs
l = sorted(zip(toas, freq, res, err), key=lambda x: x[0])

stoas, sfreq, sres, serr = zip(*sorted(zip(toas, freq, res, err)))

# Locate start and end days for each epoch
rtoas=np.around(stoas,1)
ftoas=np.floor(stoas)
ctoas=np.ceil(stoas)

epochloc=[]
for daypair in dmxranges:
    startloc=positionFirst(np.floor(daypair[0]),ftoas)
    endloc=positionLast(np.ceil(daypair[1]),ctoas)
    locpair=[startloc,endloc]
    epochloc.append(locpair)

faultypair=None
noneloc=None
for locpair in epochloc:
	if None in locpair:
		faultypair=locpair
		noneloc=locpair.index(None)

if faultypair != None:
	faultyloc=epochloc.index(faultypair)
if noneloc == 0:
	epochloc[faultyloc][0]=epochloc[faultyloc-1][1]
elif noneloc == 1:
	epochloc[faultyloc][1]=epochloc[faultyloc+1][0]
	
# Create list of 1/freq^2
invfreq=[]
for freq in sfreq:
    new=1.0/freq**2
    invfreq.append(new)

# Divide quantities by epoch
freqepochs=epochDivider(epochloc,sfreq)
invfreqepochs=epochDivider(epochloc,invfreq)
resepochs=epochDivider(epochloc,sres)
errepochs=epochDivider(epochloc,serr)

if 'crossband' in modes or 'Crossband' in modes:
	# Generate crossband points
	mathinputdir=setMode('crossband')
	mode='crossband'

	# Write quantities to files
	data=[invfreqepochs,resepochs,errepochs]
	
	n = 0
	for y in range(0,len(invfreqepochs)):
		if n < 10:
			with open(mathinputdir+"/epoch" + str(0) + str(n) +".txt", "w") as f:
				varlist = []
				for x in range(0,3):
					var = data[x][y]
					varlist.append(var)
				for i in zip(*varlist):
					f.write("{0}\t{1}\t{2}\n".format(*i))
			f.close()
			n=n+1
		else:
			with open(mathinputdir+"/epoch" + str(n) +".txt", "w") as f:
				varlist = []
				for x in range(0,3):
					var = data[x][y]
					varlist.append(var)
				for i in zip(*varlist):
					f.write("{0}\t{1}\t{2}\n".format(*i))
			f.close()
			n=n+1
	
	path='/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Output/'+fullpulsar
	if not os.path.exists(path):
		os.makedirs('/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Output/'+fullpulsar)
	
	# Write Mathematica script
	line1='files' + '= FileNames["*.txt",'+'"'+mathinputdir+'"'+'];'
	line2='data = Import[#, "Data"] & /@ files;'
	line3='m = OpenWrite["/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Output/'+fullpulsar+'/'+mode+'.txt"];'
	line4='Do[epoch = data[[i]];'
	line5='freq = epoch[[All, 1]];'
	line6='res = epoch[[All, 2]];'
	line7='err = epoch[[All, 3]]*10^-6;'
	line8='data1 = Transpose@{freq, res};'
	line9='fit = LinearModelFit[data1, x, x, Weights -> 1/err^2,VarianceEstimatorFunction -> (1 &)];'
	line10='fit2 = Normal[LinearModelFit[data1, x, x, Weights -> 1/err^2,VarianceEstimatorFunction -> (1 &)]];'
	line11='f = 0.000241022;'
	line12='t = Expand[f*fit2];'
	line13='error = f*fit["ParameterErrors"];'
	line14='error2 = Last[error];'
	line15='Write[m, {D[t, x]*-1, error2}], {i,'+str(len(epochdays))+'}];'
	line16='Close[m]'
	
	lines=[]
	for n in range(1,17):
		line=eval('line'+str(n))
		lines.append(line)
		
	with open(mathinputdir+"/temp.txt", "w") as f:
		for item in lines:
			f.write("%s\n" % item)
	
	# Run Mathematica script   
	os.chdir(mathinputdir)
	call(['math','-script','temp.txt'])

if 'dualband' in modes or 'Dualband' in modes or 'inband' in modes or 'Inband' in modes: # Need this for dualband and inband
	# Determine location of frequency bands
	freqloc=[]
	for epoch in freqepochs:
		f1=[]
		f2=[]
		f3=[]
		for n in range(len(epoch)):
			f=epoch[n]
			if obsfreqs==['820 MHz','1500 MHz']:
				if f < 1000:
					f1.append(n)
				else:
					f2.append(n)
			if obsfreqs==['430 MHz','1410 MHz']:
				if f < 800:
					f1.append(n)
				else:
					f2.append(n)
			if obsfreqs==['1410 MHz','2030 MHz']:
				if f < 1761:
					f1.append(n)
				else:
					f2.append(n)
			if len(obsfreqs)==3:
				if f < 800:
					f1.append(n)
				elif f < 1761:
					f2.append(n)
				else:
					f3.append(n)
		if len(obsfreqs)==3:
			epochfreqloc=[f1,f2,f3]
		else:
			epochfreqloc=[f1,f2]
		freqloc.append(epochfreqloc)
	
	# Arrange data based on frequency bands
	data=[freqepochs,resepochs,errepochs]
	newdata=[]
	for d in data:
		newdataepochs=[]
		for n in range(len(d)):
			epoch=d[n]
			locepoch=freqloc[n]
			newepoch=[]
			for i in locepoch:
				flist=[]
				for s in i:
					f=epoch[s]
					flist.append(f)
				newepoch.append(flist)
			newdataepochs.append(newepoch)
		newdata.append(newdataepochs)

if 'dualband' in modes or 'Dualband' in modes:
	# Generate dual band/triple band points
	mathinputdir=setMode('dualband')
	mode='dualband'
	
	# Average data
	dataav=[]
	for d in newdata:
		dav=[]
		for epoch in d:
			epochav=[]
			for values in epoch:
				av=np.mean(values)
				epochav.append(av)
			dav.append(epochav)
		dataav.append(dav)
	
	# Calculate 1/freq^2
	freqav=dataav[0]
	invfreqav=[]
	for epoch in freqav:
		newepoch=[]
		for freq in epoch:
			invfreq=1.0/freq**2
			newepoch.append(invfreq)
		invfreqav.append(newepoch)
	dataav[0]=invfreqav
	
	# Create columns
	columns=[]
	for x in range(3):
		for z in range(len(obsfreqs)):
			column=[]
			for y in range(len(epochdays)):
				value=dataav[x][y][z]
				column.append(value)
			columns.append(column)
	
	# Write columns to file
	os.chdir("/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Input/"+pulsar+"/dualband/")
	with open("dualband.csv", "wb") as f:
		writer = csv.writer(f)
		writer.writerows(columns)
	f.close()
	
	# Write Mathematica script for 3 band pulsar
	if len(epochfreqloc)==3:
		line1= 'data = Import['+'"'+mathinputdir+'/dualband.csv", "Data"];'
		line2= 's = OpenWrite["/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Output/'+fullpulsar+'/'+mode+'.txt"];'
		line3= 'f1avlist = data[[1]];'
		line4= 'f2avlist = data[[2]];'
		line5= 'f3avlist = data[[3]];'
		line6= 'f1reslist = data[[4]];'
		line7= 'f2reslist = data[[5]];'
		line8= 'f3reslist = data[[6]];'
		line9= 'f1errlist = data[[7]]*10^-6;'
		line10= 'f2errlist = data[[8]]*10^-6;'
		line11= 'f3errlist = data[[9]]*10^-6;'
		line12= 'Do[f1av = f1avlist[[i]];'
		line13= 'f2av = f2avlist[[i]];'
		line14= 'f3av = f3avlist[[i]];'
		line15= 'f1res = f1reslist[[i]];'
		line16= 'f2res = f2reslist[[i]];'
		line17= 'f3res = f3reslist[[i]];'
		line18= 'f1err = f1errlist[[i]];'
		line19= 'f2err = f2errlist[[i]];'
		line20= 'f3err = f3errlist[[i]];'
		line21= 'k = {{{f2av, f2res}, {f3av, f3res}}, {f2err, f3err}};'
		line22= 'y = {{{f1av, f1res}, {f2av, f2res}}, {f1err, f2err}};'
		line23= 'z = {{{f1av, f1res}, {f2av, f2res}, {f3av, f3res}}, {f1err, f2err, f3err}};'
		line24= 'a = Which[f1av === "nan", k, f3av === "nan", y, f1av =!= "nan" && f3av =!= "nan", z];'
		line25= 'data2 = a[[1]]; err = a[[2]];'
		line26= 'fit = LinearModelFit[data2, x, x, Weights -> 1/err^2, VarianceEstimatorFunction -> (1 &)] // Normal;'
		line27= 'fit2 = LinearModelFit[data2, x, x, Weights -> 1/err^2, VarianceEstimatorFunction -> (1 &)] ;'
		line28= 'f = 0.000241022;'
		line29= 't = Expand[f*fit];'
		line30= 'der = D[t, x];'
		line31= 'error = Last[fit2["ParameterErrors"]];'
		line32= 'adjerr = error*f;'
		line33= 'Write[s, {der*-1, adjerr}], {i,'+str(len(epochdays))+'}];'
		line34= 'Close[s]'
	
		lines=[]
		for n in range(1,35):
			line=eval('line'+str(n))
			lines.append(line)
	
	# Write Mathematica script for 2 band pulsar
	if len(epochfreqloc)==2:
		line1= 'data = Import['+'"'+mathinputdir+'/dualband.csv", "Data"];'
		line2= 's = OpenWrite["/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Output/'+fullpulsar+'/'+mode+'.txt"];'
		line3= 'f1avlist = data[[1]];'
		line4= 'f2avlist = data[[2]];'
		line5= 'f1reslist = data[[3]];'
		line6= 'f2reslist = data[[4]];'
		line7= 'f1errlist = data[[5]]*10^-6;'
		line8= 'f2errlist = data[[6]]*10^-6;'
		line9= 'Do[f1av = f1avlist[[i]];'
		line10= 'f2av = f2avlist[[i]];'
		line11= 'f1res = f1reslist[[i]];'
		line12= 'f2res = f2reslist[[i]];'
		line13= 'f1err = f1errlist[[i]];'
		line14= 'f2err = f2errlist[[i]];'
		line15= 'data2 = {{f1av,f1res},{f2av,f2res}}; err = {f1err,f2err};'
		line16= 'fit = LinearModelFit[data2, x, x, Weights -> 1/err^2, VarianceEstimatorFunction -> (1 &)] // Normal;'
		line17= 'fit2 = LinearModelFit[data2, x, x, Weights -> 1/err^2, VarianceEstimatorFunction -> (1 &)] ;'
		line18= 'f = 0.000241022;'
		line19= 't = Expand[f*fit];'
		line20= 'der = D[t, x];'
		line21= 'error = Last[fit2["ParameterErrors"]];'
		line22= 'adjerr = error*f;'
		line23= 'Write[s, {der*-1, adjerr}], {i,'+str(len(epochdays))+'}];'
		line24= 'Close[s]'
		
		lines=[]
		for n in range(1,25):
			line=eval('line'+str(n))
			lines.append(line)
	
	with open(mathinputdir+"/temp.txt", "w") as f:
		for item in lines:
			f.write("%s\n" % item)
	
	# Run Mathematica script
	os.chdir(mathinputdir)
	call(['math','-script','temp.txt'])
	os.system('rm temp.txt')

if 'inband' in modes or 'Inband' in modes:
	# Generate inband data points
	mathinputdir=setMode('inband')
	mode='inband'
	
	# Rearrange quantities based on epoch
	flist = []
	for h in range(len(obsfreqs)):
		fnlist=[]
		for n in range(len(epochdays)):
			fepoch=[]
			for i in range(len(newdata)):
				f=newdata[i][n][h]
				fepoch.append(f)
			fnlist.append(fepoch)
		flist.append(fnlist)
	
	# Create list of available bands
	bands=[]
	for n in range(1,len(flist)+1):
		f='f'+str(n)
		bands.append(f)
	
	# Create directory for each band
	for s in range(len(bands)):
		band=bands[s]
		if not os.path.exists(mathinputdir+'/'+band):
			os.makedirs(mathinputdir+'/'+band)

	# Write data to files
		n = 0
		for x in range(0,len(epochdays)):
			if n < 10:
				with open(mathinputdir+'/'+band+"/epoch" + str(0) + str(n) +".txt", "w") as f:
					varlist = []
					for y in range(0,3):
						var = flist[s][x][y]
						varlist.append(var)
					for i in zip(*varlist):
						f.write("{0}\t{1}\t{2}\n".format(*i))
				f.close()
				n=n+1
			else:
				with open(mathinputdir+'/'+band+"/epoch" + str(n) +".txt", "w") as f:
					varlist = []
					for y in range(0,3):
						var = flist[s][x][y]
						varlist.append(var)
					for i in zip(*varlist):
						f.write("{0}\t{1}\t{2}\n".format(*i))
				f.close()
				n=n+1
	
		# Write Mathematica script
		line1= 'files = FileNames["*.txt", "'+mathinputdir+'/'+band+'"];'
		line2= 'data = Import[#, "Data"] & /@ files;'
		line3= 'm = OpenWrite["/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Output/'+fullpulsar+'/'+mode+band+'.txt"];'
		line4= 'Do[epoch = data[[i]];' 
		line5= 'freq = 1/epoch[[All, 1]]^2;' 
		line6= 'res = epoch[[All, 2]];'
		line7= 'err = epoch[[All, 3]]*10^-6;' 
		line8= 'data1 = Transpose@{freq, res};'
		line9= 'fit = LinearModelFit[data1, x, x, Weights -> 1/err^2, VarianceEstimatorFunction -> (1 &)];'
		line10= 'fit2 = Normal[LinearModelFit[data1, x, x, Weights -> 1/err^2, VarianceEstimatorFunction -> (1 &)]];'
		line11= 'f = 0.000241022;'
		line12= 't = Expand[f*fit2];'
		line13= 'error = f*fit["ParameterErrors"];'
		line14= 'error2 = Last[error];'
		line15= 'Write[m, {D[t, x]*-1, error2}], {i,'+str(len(epochdays))+'}];'
		line16= 'Close[m]'
	
		lines=[]
		for n in range(1,17):
			line=eval('line'+str(n))
			lines.append(line)
		
		with open(mathinputdir+"/temp"+band+".txt", "w") as f:
			for item in lines:
				f.write("%s\n" % item)
			
		os.chdir(mathinputdir)
		call(['math','-script','temp'+band+'.txt']) # Run Mathematica script
		os.system('rm temp'+band+' .txt') # Remove temporary file

# Create plotting lists
mathout=[]
legends=[]
colors=[]
if 'crossband' in modes or 'Crossband' in modes:
	cb='/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Output/'+fullpulsar+'/crossband.txt'
	mathout.append(cb)
	legends.append('Crossband DM')
	colors.append('blue')
if 'dualband' in modes or 'Dualband' in modes:
	db='/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Output/'+fullpulsar+'/dualband.txt'
	mathout.append(db)
	legends.append('Dualband DM')
	colors.append('red')
if 'inband' in modes or 'Inband' in modes:
	ib1='/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Output/'+fullpulsar+'/inbandf1.txt'
	mathout.append(ib1)
	legends.append('Inband DM #1 '+ obsfreqs[0])
	colors.append('green')
	ib2='/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Output/'+fullpulsar+'/inbandf2.txt'
	mathout.append(ib2)
	legends.append('Inband DM #2 '+ obsfreqs[1])
	colors.append('orange')
	if len(obsfreqs)==3:
		ib3='/home/ggrillo93/Documents/Research/NANOGrav/DM/Mathematica/Output/'+fullpulsar+'/inbandf3.txt'
		mathout.append(ib3)
		legends.append('Inband DM #3 '+ obsfreqs[2])
		colors.append('grey')

# Fix errors in Mathematica output files
newout=[]
for out in mathout:
    text=fullText(out)
    badlines=[]
    for line in text:
        if not line.startswith('{'):
            badlines.append(line)
    for bline in badlines:
        text.remove(bline)
    for n in range(len(text)):
        line=text[n]
        if line.startswith('{-0.000241022*('):
            text[n]='{1000, 0}\n'
    newout.append(text)

olddm=[]
for out in newout:
    new=columnToList2(out,0)
    olddm.append(new)

newdm=[]
for dmlist in olddm:
    new=strToFloat(formatted(dmlist))
    newdm.append(new)

olddmerr=[]
for out in newout:
    new=columnToList2(out,1)
    olddmerr.append(new)
        
newdmerr=[]
for dmerr in olddmerr:
    new=strToFloat(formatted(formatted2(dmerr)))
    newdmerr.append(new)

dmxerr=np.multiply(strToFloat(columnToList2(dmxfile,2)),1000)

for n in range(0,len(mathout)):
	y=np.multiply(1000,newdm[n])
	err=np.multiply(1000,newdmerr[n])
	plt.errorbar(epochdays,y,xerr=0,yerr=err,fmt='none',ecolor=colors[n])
	plt.scatter(epochdays,y,color=colors[n], label=legends[n],marker='.')

dmxp=np.multiply(1000,np.mean(dmx)-dmx)
plt.errorbar(epochdays,dmxp, xerr=0, yerr=dmxerr, fmt='none',ecolor='purple')
plt.scatter(epochdays,dmxp,color='purple',marker='.',label=r'$\overline{DMX}$'+' - DMX')
plt.legend(fontsize=10,loc=1)
plt.title(fullpulsar+ " " +r'$\Delta$'+"DM vs Date of Obs")
plt.xlabel("Day (MJD)")
plt.ylabel(r'$\Delta$'+"DM ($pc$ $cm^{-3}$)"+"x1000")
plt.xlim(epochdays[0]-50,epochdays[-1]+50)
if modes=='Crossband':
	plt.ylim(min(dmxp)-.5,max(dmxp)+.5)
else:
	plt.ylim(min(dmxp)-1,max(dmxp)+1)
plt.minorticks_on
show()