import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

import glob

data_points = [
{'R' : [56210, 106, 481.69], 'B' : [35870, 132, 262.7], 'V' : [48536, 123, 370.6]},
{'B' : [1503881, 91, 16365], 'V' : [306725, 136, 2232]},
{'R' : [404052, 181, 2196], 'B' : [206466, 231, 887.14], 'V' : [250684, 141, 1753]},
{'R' : [295190, 118, 2452], 'V' : [271791, 104, 2581]},
{'R' : [29815, 255, 109], 'V' : [33039, 303, 102.8]},
{'R' : [22904, 185, 116.7], 'V' : [31129, 274, 109.24]}
]


galaxy_flux = 1264.2
ref_mean = 2231.92

folders = ['2_ALL', '3_BV', '4_ALL', '5_RV', '6_RV', '7_RV']
for i, f in enumerate(folders):
    point = data_points[i]
    for fi in point.keys():
        file_name = "FINAL SN/%s/LIGHTCURVE%s_%s.fit" % (f, i + 2, fi)

        data = fits.open(file_name)[0].data

        mean = np.median(data)

        print(mean, ref_mean)

        # Compare this to ref_mean
        scale = mean / ref_mean
        scaled_gflux = galaxy_flux * scale

        # Subtract off from total flux
        #point[fi][0] -= scaled_gflux


print(data_points)

type_IA = """-9.888536 -14.001578
-8.174646 -14.994168
-6.632145 -15.995169
-5.089644 -17.004583
-3.204365 -17.997173
-2.518809 -18.291585
-1.147697 -18.745821
0.394804 -19.065469
1.937305 -19.141175
2.965639 -19.090704
6.907586 -18.712174
10.335366 -18.476645
15.134258 -18.299997
18.733427 -18.165409
24.560653 -17.929879
32.958714 -17.416760
41.185386 -16.937289
48.212335 -16.533523
60.723732 -15.532522
68.436237 -14.994168
82.661524 -14.001578"""

type_IA = type_IA.split('\n')
rdays  = []
rmags = []
for t in type_IA:
    s = t.split(' ')
    rdays.append(float(s[0]))
    rmags.append(float(s[1]) + 37)

dates = [11, 13, 26, 27, 28, 29]

coefficients = """T_VR = 0.857, C_VR = 0.176; T_V = 0.108, C_V = 24.79
T_BV = 0.859, C_BV = 0.52; T_B = 0.329, C_B = 27.242
T_BV = 1.00347, C_BV = -0.1463; T_B = 0.306, C_B = 24.7167
T_VR = 0.75001, C_VR = 0.1108; T_V = -0.269, C_V = 25.311
T_VR = 0.905, C_VR = 0.386; T_V = -0.217, C_V = 24.814
T_VR = 0.89588, C_VR = 0.347; T_V = -0.0976, C_V = 24.452"""
coefficients = coefficients.replace(';', ',')

split_coefficients = coefficients.split('\n')
coefficients_list = []
for s in split_coefficients:
    variable_set = {}
    s = s.split(',')
    print(s)
    for var in s:
        #print(var)
        name, val = var.split('=')
        variable_set[name.replace(' ', '')] = float(val)

    coefficients_list.append(variable_set)

mags = []

for i, p in enumerate(data_points):

  # Check if VR is present
  if coefficients_list[i].get('T_VR'):
      V, R = p['V'], p['R']
      Vflux_sub = (V[0] - V[2] * V[1])/3
      Rflux_sub = (R[0] - R[2] * R[1])/3

      # Get VR magnitudes and transform them into standard magnitude
      v, r = -2.5 * np.log10(Vflux_sub), -2.5 * np.log10(Rflux_sub)
      v_r = v - r

      clist = coefficients_list[i]

      V_R = clist['T_VR'] + clist['C_VR']
      V = v + clist['T_V'] * V_R + clist['C_V']
      R = V - V_R
      print("On day %s, V magnitude = %f" % (dates[i], V))
      mags.append(V)
      #mags_R.append(R)
  else:
      V, B = p['V'], p['B']
      Bflux_sub = (B[0] - B[2] * B[1])/3
      Vflux_sub = (V[0] - V[2] * V[1])/3
      # Get BR magnitudes and transform them into standard magnitude
      b, v = -2.5 * np.log10(Bflux_sub), -2.5 * np.log10(Vflux_sub)
      b_v = b - v

      clist = coefficients_list[i]
      B_V = clist['T_BV'] + clist['C_BV']
      B = b + clist['T_B'] * V_R + clist['C_B']

      V = B - B_V
      print("On day %s, V magnitude = %f" % (dates[i], V))
      mags.append(V)

plt.scatter(np.array(dates) - 12, np.array(mags) + 2)
plt.plot(rdays, rmags)
plt.xlabel("Days since peak brightness")
plt.ylabel("Absolute magnitude")
plt.title("Light Curve of SN2019iee")
plt.ylim(reversed(plt.ylim()))
plt.show()
