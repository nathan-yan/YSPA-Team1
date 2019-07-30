from astropy.io import fits
from astropy.table import Table

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score
files = {
    "cat"    : "3dhst/goodss_3dhst.v4.1.cat.FITS",
    "fout"   : "3dhst/goodss_3dhst.v4.1.fout.FITS",
    "zout"   : "3dhst/goodss_3dhst.v4.1.zout.FITS",
    "master" : "3dhst/goodss_3dhst.v4.1.master.RF.FITS"
}

for f in files.keys():
    files[f] = Table.read(files[f], format = 'fits').to_pandas()
    files[f].set_index('id')

filters = 'U.B.V.R.I.J.H.Ks'.lower()
filters = ['f_' + f for f in filters.split('.')]

jcat = files['cat'][filters]
jfout = files['fout'][['lage', 'lmass', 'lsfr', 'av']]
jzout = files['zout'][['z_peak', 'z_spec']]
jmaster = files['master'][['l153', 'l155', 'l161']]

def merge(tables):
    current = tables[0]

    for i in range (len(tables) - 1):
        print(current.columns)
        print(tables[i + 1].columns)
        current = current.join(tables[i + 1], how = 'outer')

    return current

master = merge([jcat, jfout, jzout, jmaster])
print(master.loc[0])

dataset = master[['l153', 'l155', 'l161', 'lsfr']]
dataset = dataset.dropna()

x = dataset[['l153', 'l155', 'l161']]
y = dataset['lsfr']

y[y > -0.5] = 1
y[y < -0.5] = 0

X, Y = x.values, y.values
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.33, random_state = 4)
enc = preprocessing.LabelEncoder()
test_scores = enc.fit_transform(Y_test)
train_scores = enc.fit_transform(Y_train)

scaler = StandardScaler()
train_feats = scaler.fit_transform(X_train)
test_feats = scaler.fit_transform(X_test)

#model = RandomForestClassifier(n_estimators = 100)  # 61.7% accuracy on test split
model = MLPClassifier()
model.fit(train_feats, train_scores)

pred = model.predict(test_feats)

print(accuracy_score(pred, test_scores))
