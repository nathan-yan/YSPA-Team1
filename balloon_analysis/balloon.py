import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
headers = "millisec, uSv/h, CPS, CPM, Temp, Alt, roll, pitch, heading, x_accel, y_accel, z_accel, balloonstate".replace(' ', '').split(',')

f = np.loadtxt('balloon-data-2018.txt', delimiter = ',')
data = {h : [] for h in headers}

for point in f:
    for v in range (len(point)):
        data[headers[v]].append(point[v])

# CLean alt
to_delete = []
alt = np.copy(data['Alt'])

millis = np.copy(data['millisec'])

for i in range (5, len(data['Alt']) - 5):
    running_avg = np.mean(alt[i - 5 : i + 5])

    if (alt[i] - running_avg) ** 2 > 100:
        to_delete.append(i)

alt = np.delete(alt, to_delete)
millis = np.delete(millis, to_delete)

plt.subplot(2, 2, 1)
plt.hist2d(millis/1000., alt, norm = LogNorm(), bins = 100)
plt.xlabel("Time (seconds)")
plt.ylabel("Altitude (meters)")
plt.title("Altitude vs Time")

temp = np.copy(data['Temp'])
millis = np.copy(data['Alt'])


indices = np.where(temp > -40)[0]

plt.subplot(2, 2, 2)
plt.hist2d(millis[indices], temp[indices], norm = LogNorm(),bins = 50)
plt.xlabel("Altitude (meters)")
plt.ylabel("Temperature (celsius)")

plt.title("Temperature vs Altitude")
rad = np.copy(data['uSv/h'])
millis = np.delete(millis, to_delete)
rad = np.delete(rad, to_delete)
indices = np.where(millis > -10000)
plt.subplot(2, 2, 3)
plt.hist2d((millis[indices])[1000:], (rad[indices])[1000:], norm = LogNorm(),bins = 50)
plt.xlabel("Altitude (meters)")
plt.ylabel("Radiation (uSv/h)")

plt.title("Radiation vs Altitude")

plt.show()
