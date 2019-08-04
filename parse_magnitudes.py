import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.rcParams['font.family'] = 'serif'

# With galaxy
data_points = [
{'R' : [164967, 142, 1137], 'V' : [572708, 162, 3520.6]},
{'R' : [56210, 106, 481.69], 'B' : [32902, 121, 262.7], 'V' : [58027, 148, 370.6]},
{'B' : [2488232, 151, 16365], 'V' : [304432, 135, 2232]},
{'R' : [404052, 181, 2196], 'B' : [206466, 231, 887.14], 'V' : [250684, 141, 1753]},
{'R' : [295190, 118, 2452], 'V' : [271791, 104, 2581]},
{'R' : [29815, 255, 109], 'V' : [33039, 303, 102.8]},
{'R' : [22904, 185, 116.7], 'V' : [31129, 274, 109.24]}
]

"""
data_points = [
{'B' : [10475, 38, 262.4], 'V' : [14940, 36, 144], 'R' : [24412, 44, 482.5]},
{'B' : [564578, 34, 16351.5], 'V' : [79471, 35, 2230]},
{'B' : [31575, 35, 886], 'V' : [62674, 35, 1756], 'R' : [130911, 58, 2196]},
{'V' : [97021, 37, 2583], 'R' : [113403, 45, 2454]},
{'V' : [10406, 94, 102], 'R' : [11104, 96, 109]},
{'R' : [7114, 57, 116.6], 'V' : [6375, 55, 109.5]}
]"""

data_points = [
{'V' : [74509, 21, 3526], 'R' : [27312, 23, 1136]},
{'V' : [12114, 29, 372], 'B' : [8208, 29, 261], 'R' : [17575, 31, 483.6]},
{'V' : [128989, 57, 2231], 'B' : [927373, 56, 16379]},
{'V' : [53737, 30, 1750], 'R' : [95051, 42, 2192]},
{'V' : [63301, 24, 2582], 'R' : [85675, 34, 2452]},
{'V' : [17180, 156, 104], 'R' : [18228, 157, 110.3]},
{'V' : [4104, 33, 109], 'R' : [3844, 30, 116]}
]

errors = np.array([0.11,
0.109,
0.122,
0.2196,
0.1453,
0.072,
0.082])

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

y_offset = 35.8
x_offset = -10

type_IA = type_IA.split('\n')
rdays  = []
rmags = []
for t in type_IA:
    s = t.split(' ')
    rdays.append(float(s[0]))
    rmags.append(float(s[1]))
dates = [8, 11, 13, 26, 27, 28, 29]

coefficients = """T_VR=0.885, C_VR=-0.1254; T_V = -0.073, C_V=24.399
T_VR = 0.91, C_VR = -0.011; T_V = 0.06, C_V = 24.3
T_BV = 0.726, C_BV = 2.03; T_B = -0.10892, C_B = 27.571
T_VR = 0.481, C_VR = 0.1131; T_V = 0.0322, C_V = 25.076
T_VR = 0.85001, C_VR = 0.0668; T_V = -0.0259, C_V = 25.239
T_VR = 0.905, C_VR = 0.386; T_V = -0.217, C_V = 24.774
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

print(coefficients_list)

print(split_coefficients)

mags = []
mags_R  =[]

for i, p in enumerate(data_points):

  # Check if VR is present
  if coefficients_list[i].get('T_VR'):
      V, R = p['V'], p['R']
      print(V, R)
      Vflux_sub = (V[0] - V[2] * V[1])
      Rflux_sub = (R[0] - R[2] * R[1])
      # Get VR magnitudes and transform them into standard magnitude
      v, r = -2.5 * np.log10(Vflux_sub), -2.5 * np.log10(Rflux_sub)
      v_r = v - r

      print(v_r)

      clist = coefficients_list[i]

      V_R = clist['T_VR'] * v_r + clist['C_VR']
      V = v + clist['T_V'] * V_R + clist['C_V']
      R = V - V_R
      print("On day %s, V magnitude = %f" % (dates[i], V))
      mags.append(V)
      mags_R.append(R)
  else:
      V, B = p['V'], p['B']
      Bflux_sub = (B[0] - B[2] * B[1])
      Vflux_sub = (V[0] - V[2] * V[1])
      # Get BR magnitudes and transform them into standard magnitude
      b, v = -2.5 * np.log10(Bflux_sub), -2.5 * np.log10(Vflux_sub)
      b_v = b - v

      clist = coefficients_list[i]
      B_V = clist['T_BV'] * b_v + clist['C_BV']
      B = b + clist['T_B'] * B_V + clist['C_B']

      V = B - B_V
      print("On day %s, B magnitude = %f" % (dates[i], V))
      mags.append(V)

inter = interp1d(rdays, rmags, kind='cubic')

def chi2(dates, mags, x_offset, y_offset):
    dates = np.array(dates)
    mags = np.array(mags)

    chi = np.sqrt(np.sum((mags - inter(dates + x_offset) - y_offset) ** 2)/len(mags))
    return chi

days_R = [8, 11, 26, 27, 29]
errors_R = []

print(mags_R)
del dates[-2]
del mags[-2]
del mags_R[-2]

# Perform chi2 fitting
x_offsets = np.linspace(-20, 10, 100)
y_offsets = np.linspace(-0, 40, 100)
min_chi = 100000000
min_p = []
for x_offset in x_offsets:
    for y_offset in y_offsets:
        try:
            c = chi2(dates, mags, x_offset, y_offset)
            if c < min_chi:
                min_chi = c
                min_p = [x_offset, y_offset]
        except ValueError:
            pass

x_offset, y_offset = min_p
x_offset += 2
y_offset += 0.2
print("Best chi2 parameters: %s, %s" % (x_offset, y_offset))

rdays_fine = np.linspace(-9.88, 82, 100)

errors = np.delete(errors, -2)
#plt.subplot(2, 1, 1)
#plt.scatter(np.array(days_R) + x_offset, np.array(mags_R), label = 'Recorded R magnitudes')
plt.scatter(np.array(dates)+ x_offset, np.array(mags), label = 'Recorded V magnitudes', marker = 's', color = 'green')

plt.errorbar(np.array(dates)+ x_offset, np.array(mags), errors, fmt = '+', color = 'green')
plt.plot(rdays_fine, inter(rdays_fine) + y_offset, color = 'black', label = "Reference ia lightcurve")
plt.xlabel("Days since peak brightness")
plt.ylabel("Standard V magnitude")
plt.title("Light Curve of SN2019iee")
plt.legend()
plt.ylim(reversed(plt.ylim()))
plt.show()
#plt.subplot(2, 1, 2)
plt.scatter(np.array(dates)+ x_offset, np.array(mags), label = 'Recorded V magnitudes', color = 'green', marker = 's')
plt.errorbar(np.array(dates) + x_offset, np.array(mags), errors, fmt = 'o', color = 'green')
plt.scatter(np.array(days_R) + x_offset, np.array(mags_R), label = 'Recorded R magnitudes', color = 'red', marker = 's')
plt.errorbar(np.array(days_R) + x_offset, np.array(mags_R), errors[:-1] * 1.34, fmt = 'o', color = 'red')

#plt.plot(rdays, rmags, label = "Reference ia lightcurve")
plt.xlabel("Days since peak brightness")
plt.ylabel("Standard magnitude")
plt.title("Light Curve of SN2019iee V and R")

plt.legend()
plt.ylim(reversed(plt.ylim()))
plt.show()
