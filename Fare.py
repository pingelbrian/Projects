import numpy as np
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import random
import warnings
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_val_score
from sklearn import linear_model

warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")

f = open('trip_data/trip_data_1.csv','r')
g = open('trip_fare/trip_fare_1.csv','r')

F = f.readline()
G = g.readline()

setSize = 14776615

L = linear_model.RidgeCV(fit_intercept = True, cv = 4)
X = np.zeros(shape = (setSize,10))
Y = np.zeros(shape = (setSize,1))
#Load in the desired data from the dataset
for i in range(setSize):
	F = f.readline()
	G = g.readline()
	M = F.split(",")
	N = G.split(",")
	#The DateTime needs to be split and times converted into minutes from 00:00 and days of the week are binary values
	DateTime = M[5].split(' ')
	Date = DateTime[0]
	Time = DateTime[1]
	SplitDate = Date.split('-')
	DayoWeek = datetime.date(int(SplitDate[0]), int(SplitDate[1]), int(SplitDate[2])).weekday()
	TimeSplit = Time.split(':')
	Time = (int(TimeSplit[0]) * 60) + int(TimeSplit[1])
	X[i,0] = Time
	X[i,1] = M[10]
	X[i,2] = M[11]
	Y[i,:] = float(N[5])
	X[i,3 + DayoWeek] = 1
#Sample from a set size of values from the main dataset
SampleX = np.zeros(shape = (10000,10))
SampleY = np.zeros(shape = (10000,))
RandX = np.random.choice(setSize-1, size = 10000, replace = False)
print X[33,:]
for m in range(10000):
	SampleX[m,:] = X[RandX[m],:]
	SampleY[m] = Y[RandX[m],:]
Sum = 0
for n in range(10000):
	Sum = Sum + SampleY[n]

#Figure out the average fare price for comparison

Sum = Sum / 10000
print Sum
L.fit(SampleX,SampleY)


np.set_printoptions(threshold = np.nan)
#Set boundaries of the coordinate system
Increment = 300

North = 40.917577
South = 40.477399
West = -74.25909
East = -73.700009

#Find the distance that we'll generate test points over and define an increment

WtoEInc = abs(West - East)
WtoEInc = WtoEInc/Increment
NtoSInc = abs(North - South)
NtoSInc = NtoSInc/Increment
CoordX = np.zeros(shape = (Increment, 1))
CoordY = np.zeros(shape = (Increment, 1))
Locations = np.zeros(shape = (Increment * Increment, 2))

#initiate a coordinates pair matrix to feed into the model
for Lat in range(Increment):
	CoordY[Lat] = North - (Lat * NtoSInc)
	for Long in range(Increment):
		Locations[Lat * Increment + Long, 0] = West + (Long * WtoEInc)
		Locations[Lat * Increment + Long, 1] = North - (Lat * NtoSInc)
for X in range(Increment):
	CoordX[X] = West + (X * WtoEInc)
Tester = np.zeros(shape = (Increment * Increment, 10))
#720
#900
#Create a set of data to feed into the model and make predictions on

#Set time as minutes after midnight, time is variable to test times of day
Tester[:,0] = 720
Tester[:,1] = Locations[:,0]
Tester[:,2] = Locations[:,1]
Tester[:,5] = 1
prediction1 = L.predict(Tester)
fixed1 = np.zeros(shape = (Increment * Increment,))
z_min = prediction1.min()
z_max = prediction1.max()
C = np.reshape(prediction1,(Increment, Increment))
xv, yv = np.meshgrid(CoordX,CoordY)
plot_name = "heatmap 1"


color_map = plt.cm.gist_heat #plt.cm.rainbow #plt.cm.hot #plt.cm.gist_heat
ax = sns.heatmap(C, vmin = z_min, vmax = z_max, cmap = color_map)
plt.title(plot_name)
plt.show()

