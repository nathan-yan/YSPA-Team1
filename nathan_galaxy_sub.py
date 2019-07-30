import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits

import astropy
from astropy import units as u
from astropy.coordinates import Angle
from astropy.table import Table
from astrophysics_primitives import image
from astropy.modeling import models,fitting
import warnings


data = image.Image("SN2019iee_data/july12/comb/SN2019IEE_V_comb.new", y_flip = False)
data.data -= 2240

# Get pixel coordinates of the supernova
supernova_coords = [[17 * 15 + 5./60 * 15 + 10.2/3600. * 15, 41 + 46/60. + 46/3600.]]

xy = data.get_pix(supernova_coords)[0]

galaxy_cutout = np.ma.array(data.data[499:507, 500:510])
galaxy_cutout[1:4, 6:9] = np.ma.masked



def func(x, a, c, d):
    return a*np.exp(c*x)+d

means = np.mean(galaxy_cutout, axis = 0)

curve_cut = means[:5]
from scipy.optimize import curve_fit
popt, pcov = curve_fit(func, np.arange(len(curve_cut)), curve_cut, p0=(400, 0.3, 2000))
x = np.arange(len(curve_cut))
best_fit_y = func(x,popt[0],popt[1],popt[2])
print(popt)

model = astropy.modeling.functional_models.Sersic2D(amplitude = 0.7, r_eff = 25, n=4, x_0=4, y_0=5,
               ellip=0.1, theta=-1)

model.ellip.fixed = True
model.theta.fixed = True
model.x_0.fixed = True
model.y_0.fixed = True

x,y = np.meshgrid(np.arange(10), np.arange(8))

print(np.median(data.data), "MEDIAN")

fit_p = fitting.LevMarLSQFitter()

with warnings.catch_warnings():
    # Ignore model linearity warning from the fitter
    warnings.simplefilter('ignore')
    p = fit_p(model, x, y, galaxy_cutout)


img = model(x, y)
img[5][4] = img[5][3] * 5.2/4.

plt.subplot(2, 2, 1)
plt.imshow(img, vmin = -40, vmax = 160)

plt.subplot(2, 2, 2)
plt.imshow(galaxy_cutout - img, vmin = -40, vmax = 160)

means = np.mean(galaxy_cutout, axis = 0)
plt.subplot(2, 2, 4)
plt.plot(means)
plt.plot()


plt.subplot(2, 2, 3)
r = np.mean(galaxy_cutout - img, axis = 0)
r = np.delete(r, 4)
means = r
plt.plot(means)

plt.show()
