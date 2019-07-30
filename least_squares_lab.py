import numpy as np
import matplotlib.pyplot as plt


def adv_least_squares(dx, dy, order):
    row = [np.sum(dx ** i) for i in range (order)]

    highest_order = order
    data_matrix = np.zeros((order, order))

    for o in range(order):
        data_matrix[o] = row[o : o + order]
        row.append(np.sum(dx ** (o + order)))

    result_vector = np.array([[np.sum(np.multiply(dx ** i, dy))] for i in range (order)])

    inverse = np.linalg.inv(data_matrix)

    coefficients = np.dot(inverse, result_vector)
    return coefficients

def linear_fit(dx, dy):
  return adv_least_squares(dx, dy, 2)

def poly_fit(dx, dy):
    return np.polyfit(dx, dy, 2)[::-1]  # np returns highest order first

def poly_fit_with_uncertainty(dx, dy):
    y_error = dy * 0.13     # assume %13 error
    coeff, cov= np.polyfit(dx, dy, 2, w = 1/y_error, cov = True)
    fit_errors = np.sqrt(np.diag(cov))


    return coeff[::-1], cov, fit_errors

def residuals_linear(dx, dy):
    b,m = linear_fit(dx, dy)

    predicted = dx * m + b

    residuals = dy - predicted
    return residuals, np.sum(residuals)

def get_outliers(data, predictions):
    sigma = std_dev(data)

    indices = np.where(abs(data - predictions) > 2 * sigma)

    return indices

if __name__ == "__main__":
    import csv

    data = {}
    with open("starcluster_data/cal_data_VB.csv", 'r') as o:
        reader = csv.reader(o)
        headers = reader.next()
        headers = [h.replace(' ', '') for h in headers]
        data = {h : [] for h in headers}

        for row in reader:
            for i, element in enumerate(row):
                data[headers[i]].append(float(element))
    print(data)
    b_flux = np.array(data['B_flux'])
    v_flux = np.array(data['V_flux'])
    B_std = np.array(data['B_std'])
    V_std = np.array(data['V_std'])

    b = -2.5 * np.log10(b_flux)
    v = -2.5 * np.log10(v_flux)

    dx = (b - v)
    dy = (B_std - V_std)

    b_, m_ = linear_fit(dx, dy)
    print("Slope: %f; Bias: %f" % (b_, m_))

    ma, mi = np.max(dx), np.min(dx)

    plt.subplot(1, 2, 1)
    plt.scatter(dx, dy)
    plt.plot([mi, ma], [b_ + mi * m_, b_ + m_ * ma])
    plt.title("Color Calibration Line")

    # Just did color calibration above
    # Now do actual magnitude transformation

    dy = B_std - b
    dx = B_std - V_std

    b_, m_ = linear_fit(dx, dy)
    print("Slope: %f; Bias: %f" % (b_, m_))

    ma, mi = np.max(dx), np.min(dx)

    plt.subplot(1, 2, 2)
    plt.scatter(dx, dy)
    plt.plot([mi, ma], [b_ + mi * m_, b_ + m_ * ma])
    plt.title("Magnitude Calibration Line")

    #plt.xlim([-1.25, 0.75])
    #plt.ylim([-0.8, 0.3])

    plt.show()

#print(poly_fit(dx, dy))
#print()
#print(poly_fit_with_uncertainty(dx, dy))
