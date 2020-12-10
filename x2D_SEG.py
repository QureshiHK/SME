###2D SEGMENTATION TOOL, FOR SERIES OF TIFF STACKS ARRANGED INTO SEPARATE FOLDERS (PER STACK). SCROLL DOWN TO __INIT__(MAIN) FOR USER INPUT PARAMETERS


import time
start=time.time()
import cv2 as cv
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import itertools
from scipy.signal import find_peaks, peak_prominences
from skimage.segmentation import flood, flood_fill
from PIL import Image, ImageDraw
import pandas as pd
from sklearn.cluster import DBSCAN
from itertools import chain
import os
import os.path
from os.path import basename
from multiprocessing import Pool, Value, Lock, Manager
from ctypes import c_char_p
from pathlib import Path
import csv
import queue

def value_check(value):
	try:
		value = float(value)
		if value%1==0:
			value=int(value)
			return value
		else:
			value=float(value)
			return value
	except:
		value = value
		return value


#def sliceBySliceparent(name,output_dir,max_nuc_diam,min_bright,peak_dist,prom,min2Darea,max2Darea,cores):
def sliceBySlice(name,output_dir,max_nuc_diam,min_bright,peak_dist,prom,min2Darea,max2Darea,cores):

	savDir = output_dir
	max_nuc_diam = value_check(max_nuc_diam)
	min_bright=value_check(min_bright)
	prom = value_check(prom)
	min2Darea = value_check(min2Darea)
	max2Darea = value_check(max2Darea)
	peak_dist = value_check(peak_dist)
	cores=value_check(cores)
	'''
	savDir = output_dir
	max_nuc_diam = max_nuc_diam
	min_bright = min_bright
	prom = prom
	min2Darea = min2Darea
	max2Darea = max2Darea
	peak_dist = peak_dist
	cores = cores
	'''
	print("?!?! INPUTS ?!?!", savDir, type(savDir),
	max_nuc_diam, type(max_nuc_diam),
	min_bright, type(min_bright),
	prom, type(prom),
	min2Darea, type(min2Darea),
	max2Darea, type(max2Darea),
	peak_dist, type(peak_dist),
	cores, type(cores))
	#################
	#Detect features#
	#################
	cur_path=str(os.getcwd())
	parent_path=os.path.dirname(cur_path)
	savDir=output_dir+"SEG_"+basename(cur_path)
	#savDirP = savDir ######################SELF#####xx
	#nameP = name  ######################SELF#####xx
	print("!!!!!!!",savDir,"!!!!!!!!!!")
	if name.endswith(".tif"):
		name = os.path.basename(name)
		print(name)
		img2 = np.array(cv.imread(name,-1))
		print("IMG2ARRAY=", img2)
		xT,yT=img2.shape
		print("xT,yT=",xT,yT)
		xx = int(xT)
		yy = int(yT)
		gausB_img2=cv.GaussianBlur(img2, (31,31),0) #denoise to make peak detection easier  gaussian blur factor as var--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------[] O []
		#img2=cv.GaussianBlur(img2, (31,31),0)
		mask=np.zeros(img2.shape,np.uint8) #(blank image has multiple uses later on)
		clahe = cv.createCLAHE(clipLimit=20,tileGridSize=(8,8)) #exaggerate contrast to facilitate watershed, keep the clip limit low, increasing it messes up the watershed to output circles . clip limit as var-------------------------------------------------------------------------------------------------------------[] O []
		img2 = clahe.apply(img2)

		nucleus_x=[]
		nucleus_y=[]
		nucleus_intense=[]
		row_len = list(range(0,xx)) # make range variable (read dimensions of image or sm)--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------[] O []
		print("before 1st for loop")
		for xRow in range(0,yy):
			row = xRow
			xVal = gausB_img2[row]
			#feature detection by way of 1D peak detection. Hone in bright spots. Noise will be picked up, but the point is to detect everything + loads of false positives and delete the features we don't want through 2 levels of environment sensing: 1 at watersheed seeding level, and 1 at post watershed level
			peaks = find_peaks(xVal,height=min_bright,width=(int(min2Darea/2),int(xx*0.75)),distance=peak_dist,prominence=prom)  #-all peak params as vars-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------[] O []
			#print("peaks = {}".format(peaks))
			peaks_x = peaks[0]
			peaks_y = peaks[1].get('peak_heights')
			#####

			prominences = peak_prominences(xVal,peaks[0])
			#print("prominences = {}".format(prominences))
			contour_heights = xVal[peaks_x] - prominences[0]
			#print(peaks_x)
			#print(len(nucleus_x[0]))
			if len(peaks_x)>0:
				for a in peaks_x:
					nucleus_x.append(a)
					nucleus_y.append(row)
					nucleus_intense.append(xVal[a])
			else:
				continue

		print("after loop")
		peak_crd = pd.DataFrame(list(zip(nucleus_x,nucleus_y,nucleus_intense)),columns=['x','y','intensity'])
		##################################
		#use detections to set up watershed seeds#
		##################################

		mask_output=np.zeros((img2.shape[0]-2,img2.shape[1]-2),np.uint8) #final post contour fill watershed segmented images are 1918x1918 by virtue of the floodfilling method. Unsegmented blanks need to but cut to size too.
		if len(peak_crd)==0: #If no peaks detected, output blank z slice and move onto next slice
			#os.chdir(output_dir_path)
			cv.imwrite(savDir+name[:-4]+"_seg.tif",mask_output)
			return
		clusterer = DBSCAN(eps=peak_dist,min_samples=5,algorithm='auto').fit(peak_crd) #cluster local peaks together.  eps and min_samples as var-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------[] O []

		peak_crd['labels'] = clusterer.labels_
		#print(peak_crd)

		def ellbound(dataframe,label,r_index):
			boundr=int(dataframe[str(label)][r_index])
			return boundr

		mid_val = peak_crd.groupby(['labels']).median()
		dbscan_labels =mid_val.index.tolist()
		#print("dbscan_labels",dbscan_labels)
		labels_values_col=list(peak_crd.columns)
		bright_window=pd.DataFrame(columns=labels_values_col)
		for labels in dbscan_labels: #Create a sliding window  algorithm  to slide through peak clusters and identify brightest region. Finding brightest single peak alone does not suffice as it may just be a noisy single peak. Find brightest window, and locate brightest point in that to represent seed for that potential feature
			labels_values=peak_crd.loc[peak_crd['labels']==labels]

			loop_timer=0

			median_list=[]
			median_iteration_list=[]
			for a in range(0,len(labels_values)-4):
				#print("loop_timer=",loop_timer)
				aa=labels_values['intensity'].iloc[a]
				ab=labels_values['intensity'].iloc[a+1]
				ac=labels_values['intensity'].iloc[a+2]
				ad=labels_values['intensity'].iloc[a+3]
				ae=labels_values['intensity'].iloc[a+4]

				#print("ae={},ab={},ac={},ad={},ae={},loop_iteration={}".format(aa,ab,ac,ad,ae,loop_timer))
				slid_win_array = np.array([aa,ab,ac,ad,ae]) #create array of  values to take median (or mean) from
				sliding_window_median = np.median(slid_win_array) #get median of current sliding window
				#print("sliding_window_median = ", sliding_window_median)

				sliding_window_max=np.amax(slid_win_array)
				median_list.append(sliding_window_median)
				median_iteration_list.append(a)

				indexes_for_new_frame = labels_values.index[labels_values['intensity']==sliding_window_max]
				#print("indexes_for_new_frame = ",indexes_for_new_frame)
				input_dict={'x':labels_values.loc[indexes_for_new_frame[0],'x'],'y':labels_values.loc[indexes_for_new_frame[0],'y'],'intensity':labels_values.loc[indexes_for_new_frame[0],'intensity'],'labels':labels_values.loc[indexes_for_new_frame[0],'labels']}


				loop_timer+=1
				if loop_timer==(len(labels_values)-4):
					#print("median_list = ", median_list)
					max_median_list=max(median_list)
					enumeration_l=[]
					list_mid_ticker=0
					for val in median_list:
						if val==max_median_list:
							enumeration_l.append(list_mid_ticker)
						list_mid_ticker+=1
					median_max_index=median_iteration_list[int(np.median(enumeration_l))]

					a = median_max_index
					aa=labels_values['intensity'].iloc[a]
					ab=labels_values['intensity'].iloc[a+1]
					ac=labels_values['intensity'].iloc[a+2]
					ad=labels_values['intensity'].iloc[a+3]
					ae=labels_values['intensity'].iloc[a+4]

					slid_win_array = np.array([aa,ab,ac,ad,ae])
					sliding_window_max=np.amax(slid_win_array)
					#print("sliding_window_max = ",sliding_window_max)
					sliding_window_df = pd.DataFrame(columns=labels_values_col)
					swd_input_dicta={'x':labels_values['x'].iloc[a],'y':labels_values['y'].iloc[a],'intensity':labels_values['intensity'].iloc[a],'labels':labels_values['labels'].iloc[a]}
					swd_input_dictb={'x':labels_values['x'].iloc[a+1],'y':labels_values['y'].iloc[a+1],'intensity':labels_values['intensity'].iloc[a+1],'labels':labels_values['labels'].iloc[a+1]}
					swd_input_dictc={'x':labels_values['x'].iloc[a+2],'y':labels_values['y'].iloc[a+2],'intensity':labels_values['intensity'].iloc[a+2],'labels':labels_values['labels'].iloc[a+2]}
					swd_input_dictd={'x':labels_values['x'].iloc[a+3],'y':labels_values['y'].iloc[a+3],'intensity':labels_values['intensity'].iloc[a+3],'labels':labels_values['labels'].iloc[a+3]}
					swd_input_dicte={'x':labels_values['x'].iloc[a+4],'y':labels_values['y'].iloc[a+4],'intensity':labels_values['intensity'].iloc[a+4],'labels':labels_values['labels'].iloc[a+4]}
					sliding_window_df=sliding_window_df.append(swd_input_dicta,ignore_index=True)
					sliding_window_df=sliding_window_df.append(swd_input_dictb,ignore_index=True)
					sliding_window_df=sliding_window_df.append(swd_input_dictc,ignore_index=True)
					sliding_window_df=sliding_window_df.append(swd_input_dictd,ignore_index=True)
					sliding_window_df=sliding_window_df.append(swd_input_dicte,ignore_index=True)

					indexes_for_new_frame = labels_values.index[labels_values['intensity']==sliding_window_max]

					input_dict={'x':labels_values.loc[indexes_for_new_frame[0],'x'],'y':labels_values.loc[indexes_for_new_frame[0],'y'],'intensity':labels_values.loc[indexes_for_new_frame[0],'intensity'],'labels':labels_values.loc[indexes_for_new_frame[0],'labels']}

					bright_window=bright_window.append(input_dict,ignore_index=True)


		mid_val = bright_window

		seed_kernel=np.ones((3,3),np.uint8)
		i = Image.new("1",(xx,yy),"black")
		wshed_seed = ImageDraw.Draw(i)
		#print(int(mid_val['y'][13]))

		#local area intensity scan
		i2 = Image.new("1",(xx,yy),"black") #set up background binary image. Using foreround detected points to generate isolation zones as 'unknown' areas.
		bkgrd_seed = ImageDraw.Draw(i2)

		wshed_seed_retain = 0
		wshed_seed_delete = 0
		for a in range(0,len(mid_val)-1): #Early environmental scan: refine watershed seeds (environment sensing in 4 positions (N, E, S, W) to determine if watershed is worth seeding.
			#print(nucleus_x[a],nucleus_y[a])
			markerdilate=3
			bk_dil=max_nuc_diam #radius of unknown exclusion zone for backgound isolation #####USER INPUT##################XX -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------[] O []

			c_y=ellbound(mid_val,'y',a)
			c_x=ellbound(mid_val,'x',a)

			wshed_seed.ellipse((c_x-markerdilate,c_y-markerdilate,c_x+markerdilate,c_y+markerdilate),'white')
			bkgrd_seed.ellipse((c_x-bk_dil,c_y-bk_dil,c_x+bk_dil,c_y+bk_dil),'white')
			#wshed_seed.point((ellbound(mid_val,'x',a),ellbound(mid_val,'y',a)),'white')

		#set up watershed seeds for watershedding
		wshed_seed = np.array(i)
		wshed_seed_th = np.array(wshed_seed)
		wshed_seed_th = np.multiply(wshed_seed_th,1).astype(np.uint8)
		#set up background exclusion zones around seeds, ready for use in creating watersheds
		bkgrd_seed = np.array(i2)
		bkgrd_seed_th = np.array(bkgrd_seed)
		bkgrd_seed_th = np.multiply(bkgrd_seed_th,1).astype(np.uint8)

		#############################################
		#watershed module - use the seeds to create watersheds#
		#############################################

		back_img2=(img2/256).astype(np.uint8)
		back_img2=cv.threshold(back_img2,1,255,cv.THRESH_BINARY)

		kernel = np.ones((3,3),np.uint8)

		back_img2=cv.dilate(back_img2[1],kernel,iterations=2)

		bkgrd_seed_th = bkgrd_seed_th*255
		back_img2=bkgrd_seed_th
		unknown=cv.subtract(back_img2,wshed_seed_th)

		_, markers = cv.connectedComponents(wshed_seed_th)
		markers = markers+10

		markers[unknown==255]=0
		eightB_img2=(img2/255).astype("uint8")

		eightB_img2=cv.cvtColor(eightB_img2,cv.COLOR_GRAY2RGB) #my images are single channel but opencv watershed requires 3 channel input

		wshed_show = cv.watershed(eightB_img2,markers)
		wshed_show_og = wshed_show

		uniq_wshed = np.unique(wshed_show)

		wshed_mask=np.zeros(img2.shape,np.uint8) #set up blank canvas for corona placement in 2nd stage of environment sensing
		corona_wshed_delta=1.1 #comparison factor between watershed and its corona -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------[] O []
		final_wshed_ct = 0 #set up counting of wsheds
		discard_wshed_ct = 0 #how many wsheds are discarded
		for a in uniq_wshed: #iterate through each watershed region
			occur=np.count_nonzero(wshed_show==a)
			#print("{} frequency = {}".format(a,occur))
			if occur > int(max2Darea) or occur < int(min2Darea): #watershed size threshold in 2D ###########USER INPUT##################XX ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------[] O []
				wshed_show = np.where(wshed_show==a,-1,wshed_show)
				pass
			wshed_mask_im=np.where(wshed_show==a,np.uint8(255),np.uint8(0))
			wshed_mask=np.where(wshed_show==a)
			#print(wshed_mask)
			wshed_intensities=img2[wshed_mask[0],wshed_mask[1]]
			#print(wshed_intensities)
			wshed_mask_dilate=cv.dilate(wshed_mask_im,kernel,iterations=10) #create corona around watershed as environmental sensor. DIlate watershed, find difference between watershed and dilated to produce corona
			wshed_corona=cv.subtract(wshed_mask_dilate,wshed_mask_im)
			wshed_corona_mask=np.where(wshed_corona>0)
			print("wshed_corona = ", wshed_corona_mask)
			wshed_corona_intensities=img2[wshed_corona_mask[0],wshed_corona_mask[1]]
			#print("wshed_corona_intensities = ", wshed_corona_intensities)
			#print("mean wshed intensity = {}, mean corona intensity = {}".format(np.mean(wshed_intensities),np.mean(wshed_corona_intensities)))
			if np.median(wshed_intensities)<(corona_wshed_delta*np.median(wshed_corona_intensities)): #compare corona intensity with watershed
				wshed_show = np.where(wshed_show==a,-1,wshed_show)
				discard_wshed_ct+=1 #count discarded wsheds ######compile into per frame csv stats file
			else:
				final_wshed_ct+=1 #if wshed kept, count it. ######compile into per frame csv stats file


		###### to draw all watersheds pre watershed removal. wshed_show_og contains all constructed watersheds before corona based removal
		eightB_img3=(img2/255).astype("uint8")
		eightB_img3=cv.cvtColor(eightB_img3,cv.COLOR_GRAY2RGB)
		wshed_show_og = wshed_show_og.astype(np.uint8)
		wshed_show_og = cv.threshold(wshed_show_og,254,255,cv.THRESH_BINARY)
		contours_func_all=cv.findContours(wshed_show_og[1],cv.RETR_LIST,cv.CHAIN_APPROX_NONE)
		contours_all=contours_func_all[0]
		contour_draw_all=cv.drawContours(eightB_img3,contours_all,-1,(255,255,255),1)
		######

		wshed_show = wshed_show.astype(np.uint8)
		wshed_show = cv.threshold(wshed_show,254,255,cv.THRESH_BINARY)
		contours_func=cv.findContours(wshed_show[1],cv.RETR_LIST,cv.CHAIN_APPROX_NONE)
		contours=contours_func[0]
		blank_canvas=np.zeros(img2.shape,np.uint8)
		contour_draw=cv.drawContours(blank_canvas,contours,-1,color=255) #draw contours on blank image

		blanc_mask=np.zeros((contour_draw.shape[0]-2,contour_draw.shape[1]-2),np.uint8)

		contour_draw_fill=cv.floodFill(blanc_mask, contour_draw, (0,0),120)
		contour_draw_fill=np.where(contour_draw_fill[1]==0,255,contour_draw_fill[1])
		contour_draw_fill=np.where(contour_draw_fill==120,0,contour_draw_fill)

		print("END OF MODULE:", savDir)
		cv.imwrite(savDir+name[:-4]+"_seg.tif",contour_draw_fill)
		'''
		contour_draw=cv.drawContours(eightB_img2,contours,-1,(255,255,255),1) #overlay contours on source
		cv.imwrite(savDir+name[:-4]+"_segline.tif",contour_draw) #output wshed_seed_th
		cv.imwrite(savDir+name[:-4]+"_wshed_seed_refined.tif",bkgrd_seed_th)
		cv.imwrite(savDir+name[:-4]+"_wshed_seed_og.tif",contour_draw_all)

		with open((savDir+name[:-4]+"feature_det.csv"), 'w') as outcsv:
			writer = csv.writer(outcsv)
			writer.writerow(["name","seeds_kept","seeds_deleted","wshed_kept","wshed_deleted"])
			writer.writerow([name[:-4],wshed_seed_retain,wshed_seed_delete,final_wshed_ct,discard_wshed_ct])
		'''
'''
def intt(self,input_dir,output_dir,max_nuc_diam,min_bright,peak_dist,prom,min2Darea,max2Darea,cores):


	####
	curdirA=input_dir
	os.chdir(curdirA)
	print("curdirA = ", curdirA)
	print("list = ", os.listdir(curdirA))
	for dir in os.listdir(curdirA):
		print(os.path.abspath(dir))
		dirb = os.path.abspath(dir)
		#print(dirb)
		print(dir)
		if os.path.isdir(dir) is True:
			print("XXY")
			output_dir = output_dir+"SEG"+dir+"/"
			print("XXY")
			print(output_dir)
			os.mkdir(output_dir)
			print("XXY")
			curdir=os.chdir(dir)
			print("XXY")
			listdirec=os.listdir(curdir)
			print("XXYZ")
			for a in listdirec:
				sliceBySlice(a)


			with Pool(int(cores)) as p:
				print("XXYZa")
				p.map(sliceBySlice,listdirec)
				print("lol")
				os.chdir(curdirA)
'''
'''
if __name__=='__main__':
	#curdirA=os.getcwd()
	#user inputs
	input_dir ="/home/danio/Documents/segg/se/"
	output_dir="/home/danio/Documents/segg/out/"
	#x= 1920
	#y= 1920
	max_nuc_diam= 50
	min_bright=100
	prom=50
	peak_dist=50
	min2Darea=50
	max2Darea=5000
	cores=4
	# end user inputs
	#segment_go = segment_2D()
	#segment_go = segment_2D_stack(input,output,max_nuc_diam,min_bright,peak_dist,prom,min2Darea,max2Darea,cores)
	#segment_go
	#intt(input,output,max_nuc_diam,min_bright,peak_dist,prom,min2Darea,max2Darea,cores)
	#segment_2D()
	#segment_go = segment_2D()
	#segment_go
	curdirA=input_dir
	os.chdir(curdirA)
	print("curdirA = ", curdirA)
	print("list = ", os.listdir(curdirA))
	for dir in os.listdir(curdirA):
		print(os.path.abspath(dir))
		dirb = os.path.abspath(dir)
		#print(dirb)
		print(dir)
		if os.path.isdir(dir) is True:
			print("XXY")
			output_dir = output_dir+"SEG"+dir+"/"
			print("XXY")
			print(output_dir)
			os.mkdir(output_dir)
			print("XXY")
			curdir=os.chdir(dir)
			print("XXY")
			listdirec=os.listdir(curdir)
			print("XXYZ")
			for a in listdirec:
				sliceBySlice(a,max_nuc_diam,min_bright,peak_dist,prom,min2Darea,max2Darea,cores)

	end = time.time()
	print("execution time = {}s".format(end-start))
'''
#gui_queue=queue.Queue()
#gui_queue.sliceBySlice(name,output_dir,max_nuc_diam,min_bright,peak_dist,prom,min2Darea,max2Darea,cores)
