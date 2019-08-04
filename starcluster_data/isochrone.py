import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
import glob
from scipy import interpolate
files = glob.glob("./Isochrones/yapsi*")
isochrones = []
for f in files[:1]:
    isochrones.append(Table.read(f, format = 'ascii'))
    print("Loading " + str(f))

isochrones = [Table.read("Isochrones/yapsi_w_X0p602357_Z0p027643.dat", format = 'ascii')]

"""
for isochrone in isochrones:
    age = 0.001
    mask = np.where(isochrone['col1'] == age)

    stars = isochrone[mask]
    plt.plot(stars['col8'], stars['col6'])

plt.ylim(reversed(plt.ylim()))
plt.show()
"""

# Load cmd data
BV = np.loadtxt("ngc225_BV")
V = np.loadtxt("ngc225_V")

print(BV, V)
"""
isochrone_data = isochrones[0]
ages  = np.unique(isochrone_data['col1'])

masks = [np.where(isochrone_data['col1'] == age) for age in ages]

for mask in masks:
    stars = isochrone_data[mask]

    plt.plot(stars['col8'], stars['col6'])

plt.scatter(BV.astype(float), V.astype(float))
plt.ylim(reversed(plt.ylim()))
"""

mask=(V<17) & (BV<2) & (-0.1<BV)#MAKE SURE TO CHOOSE VALUES APPROPRIATE FOR YOUR STAR CLUSTER
xd=BV[mask]
yd=V[mask]

for isochrone_file in files[:1]:
    isochrone_file="Isochrones/yapsi_w_X0p602357_Z0p027643.dat"
    isochrone=Table.read(isochrone_file, format='ascii') #open up one of the files

    fit=[]
    ages=np.unique(isochrone['col1'])
    for ii,age in enumerate(ages):
        single_age_isochrone=isochrone[isochrone['col1']==age]

        x=single_age_isochrone['col8']
        y=single_age_isochrone['col6']
        #get the main-sequence part of the isochrone
        for i in range(len(x)):
            try:
                if x[i+1]>x[i]:
                    break
            except IndexError:
                i=len(x)
                break
        x=x[:i]
        y=y[:i]
        yoff=np.mean(y)-np.nanmean(yd)
        g=interpolate.interp1d(x,y)
        mask=((xd<=np.max(x)) & (xd>=np.min(x)))
        xdm=xd[mask]
        ydm=yd[mask]+yoff
        if len(ydm)==0:
            chi2=np.inf
            minOffset = 0
        else:
            # Change the offset
            yoffsets = np.linspace(-30, 20, 50)
            minchi2 = 10000000000
            minOffset = 0
            for yoffset in yoffsets:
                chi2=np.sum((np.clip(ydm-(g(xdm) + yoffset), -10, 10))**2)/(len(ydm)-1)
                if chi2 < minchi2:
                    minchi2 = chi2
                    minOffset = yoffset

            chi2 = minchi2

        #plt.scatter(xdm, ydm)
        #plt.scatter(xdm, g(xdm) + minOffset)
        #plt.ylim(reversed(plt.ylim()))

        #print(minOffset)

        #plt.show()

        fit.append([isochrone_file,age,chi2, minOffset - yoff])
        #print("%.2f%% finished"%(100.0*(ii+1)/len(ages)))

    ichi2min=np.argmin(np.array(fit)[:,2])
    print('best fit file name | age | chi2/(N-1)')
    print(fit[ichi2min])

    offset = fit[ichi2min][-1]
    d = 10 ** ((offset + 5)/5.)
    print("distance = %s" % d)

mask = np.where(isochrone['col1'] == fit[ichi2min][1])
data = isochrone[mask]

plt.scatter(BV, V)
plt.plot(data['col8'], data['col6'] + fit[ichi2min][-1], 'r')
plt.ylim(reversed(plt.ylim()))
plt.title("CMD of NGC225 and Best-Fit Isochrone")
plt.xlabel('Standard B - V')
plt.ylabel('Standard V')
plt.show()
raw_input("")
"""
isochrone = isochrones[0]

ages = np.unique(isochrone['col1'])
for age in ages:
    mask = np.where(isochrone['col1'] == 0.04)
    data = isochrone[mask]

    plt.plot(data['col8'], data['col6'] + fit[ichi2min][-1])

plt.scatter(xd, yd)
plt.ylim(reversed(plt.ylim()))
plt.show()
"""
