#script to overlay tracks with a MAX PROJECTION IN TIME (TRUE 4D IMAGE) to validate track production.

import mayavi.mlab as mlab
from scipy import array, misc
import numpy as np
import os
import cv2 as cv
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image, ImageDraw
import ast


def test_plot3d():

	directory = "/mnt/c/Users/hkqur/Desktop/4DmaxProj/4dmaxproj/"
	#data = np.array([plt.imread(os.path.join(directory, f)) for f in os.listdir(directory)])
	
	data = [] #load in 4D max projection - projected in time rather than in z, tiff stack (slice by slice)
	for f in os.listdir(directory):
		print(f)
		if f.endswith(".tif"):
			print("A")
			data.append(np.array(plt.imread(str(directory)+f)))
		else:
			pass
		#print(data)
	#print((data[1]))
	#print(len(data[0]))
	#print(data[0].shape)
	#data = data[:319]
	mitosisData = pd.read_csv("tracks_.csv") #load in track output  file from SME tracking pipeline. ensure this file is in the 4D max proj data directory

	mitosisSUB=mitosisData[mitosisData["GEN"].str.match("\[0")] #take your generation of interest. reduce dataframe to contain only nuclei of yoru genreation of interest, and those that are registered as dividing (i.e the ones which have a generation + 1 encoded in their track.
	
	#mitosisSUB=mitosisSUB[mitosisSUB['GEN'].str.contains('2')] #find if it is tracked to mitosis. this is just going to be 
	
	#generation of interest+1 at the end of the string/list.  these 2 values are important as they are modifiable parameters which decide which generation you are finding the time point for mitosos for.

	mitTimeList=[]
	for a in range(0,len(mitosisSUB)):
		valueInList=mitosisSUB["GEN"].iloc[a]
		genToList=ast.literal_eval(valueInList)
		xLOC=mitosisSUB["X"].iloc[a]
		xVAL=ast.literal_eval(xLOC)
	#	xVAL[:]=[x/3 for x in xVAL]
		yLOC=mitosisSUB["Y"].iloc[a]
		yVAL=ast.literal_eval(yLOC)
	#	yVAL[:]=[y/3 for y in yVAL]
		zLOC=mitosisSUB["Z"].iloc[a]
		zVAL=ast.literal_eval(zLOC)
	#	zVAL[:]=[z/3 for z in zVAL] #there seems to be no need to rescale the full size tracking data with the 1/3 size 3D image. The plot seems to try to scale them with each other, although it may not quite be perfect. The output gives an approximation for validation.
		tLOC=mitosisSUB["T"].iloc[a]
		tVAL=ast.literal_eval(tLOC)
		genLEN=len(genToList)
		#print("{},{},{},{}".format(xVAL[-1],yVAL[-1],zVAL[-1],tVAL[-1],genLEN))	

		s1 = [1]*len(xVAL)
		
		
		l = mlab.plot3d(zVAL, xVAL, yVAL, s1, tube_radius=0.5, colormap='Spectral') #the x,y,z dimensions are scrambled, between the plots. be mindful of this. For the tracks to overlay with the stack, xyz needed to be inverted.
	l = mlab.pipeline.volume(mlab.pipeline.scalar_field(data),vmin=450,vmax=800)
	#l=mlab.pipeline.scalar_field(data)
	#mlab.axes.label_format='%.2f'
	mlab.axes(xlabel='z',ylabel='y',zlabel='x',ranges=(0,1000,0,1000,0,1000),nb_labels=5) #set axes appropriate for your data dimensions
	#l = mlab.pipeline.volume((data))	
		
	return l
	
test_plot3d()
mlab.show()