from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry, CircularAperture, CircularAnnulus, DAOStarFinder
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

def stats(img):
    img = img.flatten()
    return np.std(img), np.median(img)

def color(bv):
    return  2.8E-3 / (4600 * (1./(0.92 * bv + 1.7) + 1./(0.92 * bv + 0.62)))

"""
with fits.open("M36_B.fits") as HDU:
  dataB = HDU[0].data

with fits.open("M36_V.fits") as HDU:
    dataV = HDU[0].data

with fits.open("m36-R.fit") as HDU:
    dataR = HDU[0].data
"""


with fits.open("starcluster_data/ngc225_B_comb.fit") as HDU:
  dataB = HDU[0].data

with fits.open("starcluster_data/ngc225_V_comb.fit") as HDU:
    dataV = HDU[0].data

with fits.open("starcluster_data/ngc225_V_comb.fit") as HDU:
    dataR = HDU[0].data


data = [dataB, dataV, dataR]
_, medians, sigmas = np.array([sigma_clipped_stats(d, sigma = 2.0) for d in data]).T
starFinders = [DAOStarFinder(threshold = 3.0 * sigmas[i], fwhm = 4.0) for i in range (len(data))]

stars = [starFinders[i](data[i] - medians[i]) for i in range (len(data))]
#print(stars[0])
starsBCentroids = np.array([stars[0]['xcentroid'], stars[0]['ycentroid']]).T
starsVCentroids = np.array([stars[1]['xcentroid'], stars[1]['ycentroid']]).T
print(len(starsBCentroids))
d = np.sqrt(10)

starCoordinates = []

for centroid in starsBCentroids:
    for other in starsVCentroids:
        if np.sqrt(np.sum((centroid - other)**2)) < 3:
            starCoordinates.append((other + centroid)/2.)
            break;
print(len(starCoordinates))
starCoordinates=np.array(starCoordinates)

apertureAp = CircularAperture([starCoordinates[:, 0], starCoordinates[:, 1]], r = 8)
apertureAnn = CircularAnnulus([starCoordinates[:, 0], starCoordinates[:, 1]], r_in = 9, r_out = 11)

area = apertureAnn.area()

phot_tablesB = {'aperture_sum_0': aperture_photometry(dataB - medians[0], [apertureAp]), 'aperture_sum_1':aperture_photometry(dataB - medians[0], [apertureAnn])}
phot_tablesV = {'aperture_sum_0': aperture_photometry(dataV - medians[1], [apertureAp]), 'aperture_sum_1':aperture_photometry(dataV - medians[1], [apertureAnn])}

background_sums = np.array([phot_tablesB['aperture_sum_1']['aperture_sum'], phot_tablesV['aperture_sum_1']['aperture_sum']]) / float(area)

fluxesB = np.array([phot_tablesB['aperture_sum_0']['aperture_sum'], phot_tablesB['aperture_sum_0']['aperture_sum']])
fluxesV = np.array([phot_tablesV['aperture_sum_0']['aperture_sum'], phot_tablesV['aperture_sum_0']['aperture_sum']])

fluxesB[1] -= background_sums[0] * apertureAp.area()
fluxesV[1] -= background_sums[1] * apertureAp.area()

print(fluxesV[1])

magB = -2.5 * np.log10(fluxesB)
magV = -2.5 * np.log10(fluxesV)
print(len(magV[0]), len(magV[1]))

#plt.imshow(dataB, vmin = 400, vmax = 1800)
#apertureAnn.plot()
#apertureAp.plot()

"""
for i in range (len(starCoordinates)):
    x, y = starCoordinates[i]
    print(magV[1][i])
    try:
        plt.text(x, y, str(int(fluxesB[1][i])) + " " + str(int(fluxesV[1][i])))

    except:
        pass
"""
#plt.ylim(reversed(plt.ylim()))
#plt.xlim(reversed(plt.xlim()))
plt.show()

zams_bv = [-0.3, -0.2, -0.15,-0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 , 0.7, 0.8, 0.9, 1.0, 1.2, 1.3]
zams_mag = [-3.50, -1.30, -0.50, 0.30, 1.30, 1.80, 2.25, 2.8, 3.35, 4.05, 4.60, 5.20, 5.70, 6.10, 6.60, 7.45, 7.90]

plt.subplot(1, 2, 1)
plt.plot(magB[0] - magV[0], magV[0], 'ro')

plt.xlabel("B - V")
plt.ylabel("Apparent V")
plt.title("CMD without Annulus")

#plt.xlim(-0.1, 1.7)

plt.ylim(reversed(plt.ylim()))

bv = magB[1] - magV[1]

plt.subplot(1, 2, 2)
#plt.plot(bv, magV[1] + 22.98, 'o')

BV = bv * 0.865 - 0.161

B = magB[1] + (-0.1356 * BV) + 22.98
V = B - BV

plt.plot(BV, V, 'co')

plt.xlabel("Standard B - V")
plt.ylabel("Standard V")
plt.title("CMD with Annulus")

plt.xlim(-0.25, 2.)
plt.ylim(reversed([8., 17.5]))

print("Saving CMD data")
np.savetxt("ngc225_BV", BV)
np.savetxt("ngc225_V", V)

plt.show()

plt.imshow(dataB - medians[0])
#apertureAp.plot()
#apertureAnn.plot()
plt.show()
