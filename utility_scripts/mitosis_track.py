#script to track where mitosis occurs per generations. Construct 3D graph (in 2 dimensions, 3rd dimension depicted as pointcolour. plot position of mitosis events in 2D space, colour code with time. Using shades of grey for ease of viewing.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast
from matplotlib import cm
import math

tLENHIGH=115 #upper limit of cell cycle lengths to plot (mitosis points from nuclei greater than this cell cycle length (in raw frames) are excluded. Decided by predetermined exclusion criteria)
plotxAx=0
tLENLOW=35 #lower limit of cell cycle lengths to plot (mitosis points from nuclei below this cell cycle length (in raw frames) are excluded. Decided by predetermined exclusion criteria)
plotxAx=0 #Axes to plot in 2D. plotxAx is the X axis of the scatter plot plane. The number denotes which axis of the embryo to plot. x=0 y=1 z=2
plotyAx=1

mitosisData = pd.read_csv("tracks_.csv")

mitosisSUB=mitosisData[mitosisData["GEN"].str.match("\[1")] #take your generation of interest. reduce dataframe to contain only nuclei of yoru genreation of interest, and those that are registered as dividing (i.e the ones which have a generation + 1 encoded in their track.
mitosisSUB=mitosisSUB[mitosisSUB['GEN'].str.contains('2')] #find if it is tracked to mitosis. this is just going to be generation of interest+1 at the end of the string/list.  these 2 values are important as they are modifiable parameters which decide which generation you are finding the time point for mitosos for.

mitTimeList=[]
for a in range(0,len(mitosisSUB)):
	valueInList=mitosisSUB["GEN"].iloc[a]
	genToList=ast.literal_eval(valueInList)
	xLOC=mitosisSUB["X"].iloc[a]
	xVAL=ast.literal_eval(xLOC)
	yLOC=mitosisSUB["Y"].iloc[a]
	yVAL=ast.literal_eval(yLOC)
	zLOC=mitosisSUB["Z"].iloc[a]
	zVAL=ast.literal_eval(zLOC)
	tLOC=mitosisSUB["T"].iloc[a]
	tVAL=ast.literal_eval(tLOC)
	genLEN=len(genToList)
	print("{},{},{},{}".format(xVAL[-1],yVAL[-1],zVAL[-1],tVAL[-1],genLEN))
	tLEN=tVAL[-1]-tVAL[1]+1
	
	plotcrd=(xVAL[-1],yVAL[-1],zVAL[-1],int(tVAL[-1]),genLEN,tLEN) #create tuple with all the coordinates we could need, parameters of plotcrd tuple are: x,y,z,T,gen,cellcyclelength
	mitTimeList.append(plotcrd)
print(mitTimeList)

range_finder = []
for b in mitTimeList:
	if b[5] > tLENLOW or b[5] < tLENHIGH:
		range_finder.append(b[3])
	else:
		continue
rangemin=min(range_finder)
rangemax=max(range_finder)

def rangeNorm(timevalue,rangemin,rangemax):
	#normVal=((timevalue-(rangemin-1))/(rangemax-(rangemin-1)))/254
	normVal=((timevalue-rangemin)/(rangemax-rangemin))
	print("normVal",normVal)
	return normVal
	

for b in mitTimeList:
	if b[5] < tLENLOW or b[5] > tLENHIGH:
		continue
	print("b= ", b)
	plt.scatter(b[plotxAx],b[plotyAx],c=cm.gist_gray(rangeNorm(b[3],rangemin,rangemax)),s=(b[5]*b[5])/10) #need to multiply the mitosis time value by 15 as the colour pallette doesn't register the small time differences. scaling up puts distance between time points and these register better on colour scale. make points nice and big to emphasise colours.
	
plt.title("Mitosis points in the transition from 256->512cell stage (XY plane, dataset: wtPii)")	
#plt.legend()
#plt.grid(True)	
#plt.colorbar()
plt.show()
