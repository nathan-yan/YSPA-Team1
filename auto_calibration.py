from astrophysics_primitives import image
from astropy.stats import sigma_clipped_stats
from astropy import units as u
from astropy.coordinates import Angle
from astropy.table import Table
from photutils import aperture_photometry, CircularAperture, CircularAnnulus, DAOStarFinder

import numpy as np

import matplotlib.pyplot as plt

from bs4 import BeautifulSoup as bs
import requests as r

plt.ion()

def apass_query(ra, dec, radius, items = 28):
    data = r.post("https://www.aavso.org/cgi-bin/apass_dr10_download.pl", data = {"ra" : Angle(ra, u.degree).to_string(unit = u.hour, sep = ':'), "dec" : Angle(dec, u.degree).to_string(unit = u.degree, sep = ':'), "radius":radius, "outtype":0})

    b = bs(data.text)

    star_data = [[] for i in range (items)]
    tds = b.find_all('td')

    for row in range (0, len(tds), items):
        if row == 0:
            header = [i.get_text() for i in tds[row : row + items]]
        else:
            points = [i.get_text() for i in tds[row : row + items]]

            for p in range(len(points)):
                val = points[p]
                if val != 'NA':
                    star_data[p].append(float(val))
                else:
                    star_data[p].append(val)

    return Table(star_data, names = header)

configuration = 'VB' # or "VR"

#field = image.Image("SN2019IEE_B_comb.new")
#images = ["SN2019IEE_V_comb.fit", "SN2019IEE_B_comb.new"]

#base = "SN2019iee_data/july 11 ALL(t21)/"
base = "starcluster_data/"
fname = "ngc225_$_comb"

field = image.Image(base + fname.replace('$', configuration[0]) + ".new")

#V, R, B

images = []
for c in configuration:
    print("Loading image: %s" % base + fname.replace('$', c) + '.fit')
    images.append(base + fname.replace('$', c) + '.fit')

assignments = {}
for f in "VBR":
    if f in configuration:
        assignments[configuration.index(f)] = f

saveas = base + "cal_data_VB.csv"

images = [image.Image(f) for f in images]

minimum, maximum = 0, 1000
done = False
while not done:
    s = raw_input("Stretch$ ")
    if s == 'q':
        done = True
    else:
        try:
            s = [int(i) for i in s.split(',')]
            minimum, maximum = s
            plt.clf()
            plt.imshow(field.data, cmap = "gray_r", vmin = minimum, vmax = maximum)
            plt.show()

        except:
            print("Invalid input")

# Get statistics
print("Getting statistics")
stats = sigma_clipped_stats(field.data, sigma = 3.0)
field.set_statistics(stats)

for i in images:
    stats = sigma_clipped_stats(i.data, sigma = 3.0)
    i.set_statistics(stats)

print("Finding stars")
num_stars = field.find_stars(10., 3.)
print("Found %d stars" % num_stars)

print("Resolving field")
middle = (field.width / 2., field.height / 2.)

center_coords = field.get_world([middle])[0]

print("Retrieving calibration stars")
cal_stars = apass_query(center_coords[0], center_coords[1], 0.3)

print("Found %d stars within 0.2 degree radius" % len(cal_stars))
print("Matching stars")

cal_coordinates = []

for cal_star in cal_stars:
    """
    try:
        plt.plot([cal_star['RA (deg)']], [cal_star['Dec (deg)']], 'o', ms = int(float(cal_star['Johnson V (V)']) ** 2 / 20.))
        plt.text(cal_star['RA (deg)'], cal_star['Dec (deg)'], cal_star['Johnson V (V)'])
    except Exception as e:
        print(e)
    """
    cal_coordinates.append([cal_star['RA (deg)'], cal_star['Dec (deg)']])

plt.xlabel("RA")
plt.ylabel("Dec")
plt.title("Calibration Stars")

plt.show()

matched_stars = field.match_stars(cal_coordinates, dist = 2)

print(matched_stars)

for star in range(len(matched_stars)):
    plt.text(matched_stars[star]['xcentroid'], matched_stars[star]['ycentroid'], str(star), color = 'red')
plt.show()

done = False
selected_stars = {}
while not done:
    s = raw_input("flux$ ")
    try:
        if s == 'q':
            done = True
        elif s[0] == 'd':
            idx = int(s.split(' ')[1])
            del selected_stars[idx]

            plt.text(star['xcentroid'], star['ycentroid'], str(idx), color = 'red')
            plt.show()
        else:

            idx = int(s)
            star = matched_stars[idx]
            selected_stars[idx] = star

            plt.text(star['xcentroid'], star['ycentroid'], str(idx), color = 'green')
            plt.show()
    except Exception as e:
        print(e)

done = False
R = 8
R_in = 10
R_out = 15
positions = [[selected_stars[s]['xcentroid'], selected_stars[s]['ycentroid']] for s in selected_stars.keys()]

while not done:
    s = raw_input("aperture$ ")
    s = s.split(' ')
    try:
        if s[0] == 'q':
            done = True
        elif s[0] == 'ap':
            R = int(s[1])
        elif s[0] == 'in':
            R_in = int(s[1])
        elif s[0] == 'out':
            R_out = int(s[1])

        apertureAp = CircularAperture(positions, r = R)
        apertureAnn = CircularAnnulus(positions, r_in = R_in, r_out = R_out)

        plt.cla()
        plt.imshow(field.data, cmap = "gray_r", vmin = minimum, vmax = maximum)
        apertureAp.plot()
        apertureAnn.plot()

        plt.show()
    except Exception as e:
        print(e)

print("Performing aperture photometry on images")
phot_tables = [aperture_photometry(f.data - f.median, [apertureAp, apertureAnn]) for f in images]
#phot_tables2 = aperture_photometry(field - field.median)

v_flux = []
b_flux = []
V_std = []
B_std = []
r_std = []
g_std = []

fluxes = []
for table in phot_tables:
    ann_mean = table['aperture_sum_1'] / apertureAnn.area()
    flux = table['aperture_sum_0'] - ann_mean * apertureAp.area()

    fluxes.append(list(flux))

# Get standard magnitudes

for i, s in enumerate(selected_stars.keys()):
    cal_star_idx = int(selected_stars[s]['idx'])
    cal_star = cal_stars[cal_star_idx]

    names = ['Johnson V (V)', 'Sloan r\' (SR)', 'Johnson B (B)', 'Sloan g\' (SG)']
    values = [cal_star[n] for n in names]
    if 'NA' not in values:
        V_std.append(cal_star[names[0]])
        r_std.append(cal_star[names[1]])
        B_std.append(cal_star[names[2]])
        g_std.append(cal_star[names[3]])

        print(cal_star)
    else:

        for f in range (len(fluxes)):
            del fluxes[f][i]

V_std, B_std, r_std, g_std = np.array(V_std).astype(float), np.array(B_std).astype(float), np.array(r_std).astype(float), np.array(g_std).astype(float)

# Convert r, B, g and V -> R
VR = 1/1.321 * (g_std - r_std) - (0.278)/1.321 * (B_std - V_std) + (0.219)/1.321
R_std = -VR + V_std

cal_data = np.stack([np.array(f) for f in fluxes] + [V_std, B_std, R_std, r_std, g_std], axis = -1)
print(cal_data, cal_data.T)
cal_data = cal_data
with open(saveas, 'w') as o:
    headers = []
    for i in range (len(images)):
        filter = assignments[i]
        headers.append(filter + "_flux")

    headers += ['V_std', 'B_std', 'R_std', 'r_std', 'g_std']
    o.write(configuration)
    o.write(",".join(headers))
    for row in cal_data:
        print(','.join(row.astype(str)))
        o.write(','.join(row.astype(str)) + '\n')

print("Saved calibration data")

plt.show()
input("")
