from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def bin_data(data, minimum, maximum, bins = 10):
    bin_edges = np.linspace(minimum, maximum, bins)

    indices = []
    for e in range (len(bin_edges) - 1):
        mi = bin_edges[e]
        ma = bin_edges[e + 1]

        indices.append(np.where((mi < data) & (data < ma)))

    return indices

sfr_full_fname = "gal_totsfr_dr7_v5_2.fits"
mass_full_fname = "totlgm_dr7_v5_2b.fit"
Z_full_fname = "gal_fiboh_dr7_v5_2.fits"

names = [sfr_full_fname, mass_full_fname, Z_full_fname]

data = [fits.open(n)[1].data for n in names]
print(len(data[0]), len(data[1]), len(data[2]))
data_indices = []

indices = np.where((data[0]['AVG'] > -99) * (data[1]['AVG'] != -1) * (data[2]['AVG'] > -99.9))[0]

data[0] = data[0][indices]
data[1] = data[1][indices]
data[2] = data[2][indices]

print(len(data[0]), len(data[1]), len(data[2]))

"""
plt.subplot(2, 2, 1)
plt.hist2d(data[1]['AVG'],data[2]['AVG'],bins=300, norm = LogNorm())
plt.colorbar()
plt.title('Mass/Metallicity relation for SDSS Galaxies')
plt.xlabel(r'log Mass [$M_\odot$]')
plt.ylabel(r'log Gas Phase Metallicities')

plt.subplot(2, 2, 2)
plt.hist2d(data[0]['AVG'],data[1]['AVG'],bins=300, norm = LogNorm())
plt.colorbar()
plt.title('SFR/Mass relation for SDSS Galaxies')
plt.xlabel(r'log SFR')
plt.ylabel(r'log Mass [$M_\odot$]')

plt.subplot(2, 2, 3)
plt.hist2d(data[0]['AVG'],data[2]['AVG'],bins=300, norm = LogNorm())
plt.colorbar()
plt.title('SFR/Metallicity relation for SDSS Galaxies')
plt.xlabel(r'log SFR')
plt.ylabel(r'log Gas Phase Metallicities')
"""

plt.ion()

sfr_check = bin_data(data[0]['AVG'], -2, 2, bins = 10)
for i, bin in enumerate(sfr_check):
    plt.subplot(2, 5, i + 1)
    metallicity = data[2]['AVG'][bin]
    mass = data[1]['AVG'][bin]
    plt.title("SFR Static, Mass vs Metallicity")
    plt.hist2d(mass, metallicity, bins = 100)

plt.show()
raw_input("")
mass_check = bin_data(data[1]['AVG'], 7, 12, bins = 10)
for i, bin in enumerate(mass_check):
    plt.subplot(2, 5, i + 1)
    sfr = data[0]['AVG'][bin]
    metallicity = data[2]['AVG'][bin]

    plt.title("Mass Static, SFR vs Metallicity")
    plt.hist2d(sfr, metallicity, bins = 100)

plt.show()
raw_input("")
metallicity_check = bin_data(data[2]['AVG'], 8, 9.5, bins = 10)
for i, bin in enumerate(metallicity_check):
    plt.subplot(2, 5, i + 1)
    mass = data[1]['AVG'][bin]
    sfr = data[0]['AVG'][bin]

    plt.title("Metallicity Static, SFR vs Mass")
    plt.hist2d(sfr, mass, bins = 100)

plt.show()
raw_input("")
