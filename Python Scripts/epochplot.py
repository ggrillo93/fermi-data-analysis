import numpy as np
from matplotlib.pyplot import *
from matplotlib import pyplot as plt
import uncertainties as unc
from uncertainties import unumpy as unp
import glob

def columnToList(text,colnumber):
    with open(text) as f:
        lista=[]
        for line in f:
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

def detectfsep(freqlist):
    fbandloc = [0]
    n = 0
    while n < len(freqlist)-1:
        f1 = freqlist[n]
        f2 = freqlist[n+1]
        d = f2-f1
        if d > 20.0:
            fbandloc.append(n+1)
        n = n+1
    fbandloc.append(-1)
    return fbandloc

def epochplotter(freq, res, err, startep, stopep, epochloc, av = False):
    # Plots frequency vs residuals for given epochs
    res = np.asarray(res)*10**6
    bnds = [epochloc[startep][0],epochloc[stopep][1]]
    ifreq, ires, ierr = freq[bnds[0]:bnds[1]+1], res[bnds[0]:bnds[1]+1], err[bnds[0]:bnds[1]+1]
    ifreq_s, ires_s, ierr_s = zip(*sorted(zip(ifreq, ires, ierr)))
    seploc = detectfsep(ifreq_s)
    bands = []
    for n in range(len(seploc)-1):
        bands.append([ifreq_s[seploc[n]:seploc[n+1]], ires_s[seploc[n]:seploc[n+1]], ierr_s[seploc[n]:seploc[n+1]]])
    if av == False:
        for n in bands:
            plt.figure()
            plt.errorbar(n[0], n[1], xerr = 0, yerr = n[2])
            #plt.scatter(n[0],n[1], color = 'blue')
            plt.ylim(min(n[1])-0.5, max(n[1])+0.5)
            plt.xlabel('Frequency (MHz)')
            plt.ylabel('Residuals (us)')
    else:
        for n in bands:
            d = {}
            for a, b in zip(n[0], n[1]):
                d.setdefault(a, []).append(b)
            avfreq, avres = [], []
            for key in d:
                avfreq.append(key)
                avres.append(sum(d[key])/len(d[key]))
            plt.figure()
            plt.scatter(avfreq,avres)
            plt.ylim(min(avres)-0.5, max(avres)+0.5)
            plt.xlabel('Frequency (MHz)')
            plt.ylabel('Mean residuals (us)')
    plt.show()
    return

data = np.loadtxt("J1713.txt")
data = data.T
freq, pre, res, toas, err = data[0], data[1], data[2], data[3], data[4]
parfile="J1713.par"
parameters=columnToList(parfile,0)
parvalue=columnToList(parfile,1)

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

epochplotter(sfreq,sres, serr,33,33,epochloc,av=False)
