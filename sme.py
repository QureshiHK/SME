import PySimpleGUIQt as pg
from os import path
import os
import csv
#import get_variable_name
import numpy as np
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
#from nucleus_track_3d_2 import *
from x2D_SEG import *
import subprocess
#import threading
import dummy_threading as threading #2020-11-03 this library allows the imported segmentation script to run. although I need to modify inputs to be int or float, as they pass into the imported saturnia package. dummy__threading is better as it executes threads sequentially. Parallelising threads in a shared memory space (via nromal threading library) causes crashes. May need to see if possible to nest threading, where individual iterations of saturnia run sequentially but multiple instances are spread across processes (parallelisation with no shared memory)
import _dummy_thread as threading_norm
import queue
#from saturnia import *
from output_transform import *
import multiprocessing.dummy as mp
import multiprocessing as mpg
import gc
from contextlib import closing
import traceback
#libraries for tracking function (track_nuc)
import time
import matplotlib.pyplot as plt ##!! dependency: ensure tkinter is installed, for matplotlib gui (apt-get install python3-tk)
matplotlib.use('Qt5Agg')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TK
import matplotlib.pylab as pylab
#import numpy
from numpy import ndarray as npd
import pandas as pd
import scipy.spatial as spatial
from scipy.spatial.distance import cdist
from sklearn import neighbors as skn
import sys







def value_check_trkparam(value):
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

#USER INPUT VARIABLES:
def track_nuc(trk_input,trk_output,frame_no, scan_rad, min_scanR, decay_rate, mit_buff):
	window.read(20)
	start = time.time()
	####USER INPUTS:
	#track_input

	data_len = int(frame_no) #length of dataset in frames/timepoints
	scan_rad_lim = value_check_trkparam(min_scanR) #the lowest scan radius to decay to.
	scan_rad_start = value_check_trkparam(scan_rad) #radius to scan within, to begin decay from
	scan_rad_decay_factor = value_check_trkparam(decay_rate) #multiplied by the cell generation and subtracted from scan_rad_start
	mitosis_buffer = int(mit_buff)
	####END USER INPUTS



	min_scan_diam = value_check_trkparam(min_scanR)
	cycles_to_plataeu = (scan_rad_start - min_scan_diam)/scan_rad_decay_factor #number of cell cycles before cell diameter doesn't reduce dramatically
	 #number of frames after mitosis that another mitosis event should not be triggered
	data_len = data_len + 1
	def narray(lst):
		lst = np.array(lst, dtype=np.float)
		return lst

	data = pd.read_csv(trk_input)
	metalist_X = data.groupby('T')['X'].apply(list)
	metalist_Y = data.groupby('T')['Y'].apply(list)
	metalist_Z = data.groupby('T')['Z'].apply(list)
	metalist_T = data.groupby('T')['T'].apply(list)

	T = data.groupby(['T']).groups.keys()
	T_len = len(data.groupby(['T']).groups.keys())

	#print("T len = ", T_len)
	T = [*T] #unpack dict_key produced earlier, as it is non indexable
	#print(T)

	TP_coord_meta = [] #list of lists - each sublist contains coords from a single time point,  the time point corresponds to the position of the sublist in this master list
	for a in T:
		TP_coord = [] #list holds coordinates from current time point (a), to later append into master list (TP_coord_meta) before this list is overwritten with data for a new frame
		while a<data_len: #for and while set up to limit frames we inspect
			print(T)
			X_array = metalist_X[a]
			Y_array = metalist_Y[a]
			Z_array = metalist_Z[a]
			T_array = metalist_T[a]

			X_array = narray(X_array)
			Y_array = narray(Y_array)
			Z_array = narray(Z_array)
			T_array = narray(T_array) #business as usual, extract time point relevant lists, and convert to numpy array
	#construct coord triple as list
			no_of_coords = len(X_array[:])
			for e in range(0,no_of_coords): #iterate through every set of coordinates, to add to time point specific list, i.e. these coordinates occur at this time point, the time point represented by its position in the meta/master list
				#### 2019/10/17: produce 6,1 array representing 4D coords and lineage identifiers, in current form: (x,y,z,T,lin,gen)
				c = []
				c.append(X_array[e])
				c.append(Y_array[e])
				c.append(Z_array[e]) #c is the coordinate sublist to add to master list: b
				c.append(T_array[e])
				c.append(000) #### append blank lineage ID (a 'family name')
				c.append(000) #### append blank generation ID (to iterate +1 with every division)
				#print("c = ", c)
				TP_coord.append(c)
			#print("first coord time point 0 =",TP_coord[:])
			break
		TP_coord_meta.append(TP_coord)
		#print("TP_coord_meta = ", TP_coord_meta[:]) #inspect master list

	def reshape_(coord):
		arroo = [0]
		arroo[0] = coord
		arroo = np.array(arroo)
		return arroo

	def closest_node(node,nodes):
		return nodes[cdist([node],nodes).argmin()] #function to find nearest coordinate to named coordinate (node) from a list of coordinates (nodes)

	###search for points within radius - return the number of points within radius, by looking at length of list
	def num_pts_within_dist(compare_pt, compare_lst, radius):
		#print("num:: compare_pt = ", compare_pt)
		#print("num:: compare_lst = ", compare_lst)
		tree = skn.KDTree(compare_lst, leaf_size = 2)
		#print("skn_tree = ", tree)
		ind = tree.query_radius(compare_pt, r=radius)
		#print("num_pts_within_dist, ind = ", ind)
		num_near_pts = len(ind[0])
		#print("num_near_pts = ", num_near_pts)
		return num_near_pts

	###search for points within radius - return the coordinates (how to factor in, re-appending of T-coordinate to this)
	def ID_pts_within_dist(compare_pt, compare_lst, radius):
		tree = skn.KDTree(compare_lst, leaf_size = 2)
		ind = tree.query_radius(compare_pt, radius, return_distance=True, sort_results=True) #### changes 2020-09-11 checking return distance
		#print("ID:: ind = ", ind)
		#print("ID:: ind[0] = ", ind[0])
		ind = ind[0]
		near_pts = []
		#print("ID:: compare_lst = ", compare_lst)
		#print("ID:: compare_pt = ", compare_pt)
		for i in ind[0]:
			#print("ID:: i = ", i)
			near_pts.append(compare_lst[i])
			#print("near_pts = ", near_pts)
		return near_pts

	def nearest_K_pts_v2(compare_pt,compare_lst, K):
		#print("nearest_K_pts(skn):: compare_pt = ", compare_pt)
		#print("nearest_K_pts(skn):: compare_lst = ", compare_lst)
		tree = skn.KDTree(compare_lst, leaf_size = 2)
		dist, ind = tree.query(compare_pt, k=K)
		#print("nearest_K_pts(skn):: ind = ", ind)
		#print("nearest_K_pts(skn):: dist = ", dist)
		near_K = []
		for i in ind[0]:
			#print("nearest_K_pts(skn):: i = ", i)
			near_K.append(compare_lst[i])
			#print("nearest_K_pts(skn):: near_K = ", near_K)
		return near_K

	def SLICE(tracked):
		track_slice = tracked[:3]
		arroo = [0]
		arroo[0] = track_slice
		return arroo

	def STICK(arroo, tracked_remain):
		#print("arroo[0] = ", arroo[0])
		if isinstance(arroo[0],list):
			print("is list!")
			f = arroo
			pass
		elif isinstance(arroo[0],float):
			print("is float!")
			f = arroo
			if isinstance(arroo,np.ndarray):
				f = npd.tolist(arroo)
			else:
				pass
			pass
		else:
			f = npd.tolist(arroo[0])
		#print("arroo = ", arroo)
		print("tracked remain = ", tracked_remain)
		f.append(tracked_remain[0])
		#print("f = ",f)
		f.append(tracked_remain[1])
		f.append(tracked_remain[2])
		#print("f =", f)
		return f

	#### THE MITOSIS ENGINE:  TRACKING USING NEAREST NEIGHBOUR PROCESS ####
	track_list = []
	#### append generation 0 to every coordinate prior to this upcoming loop, and then +1 to generation value whenever coordinate falls into a mitotic lineage.
	'''
	if  isinstance(TP_coord_meta,list):
		print("TP_coord_meta is a list")
	elif isinstance(TP_coord_meta, np.ndarray):
		print("TP_coord_meta is an array")
	else:
		print("not sure what TP_coord_meta is")
	'''
	family_ID_clock = 0
	for a_crd in TP_coord_meta[0]:          ####TP_coord_meta[0] is the list of coordinates in time point 0 but time point 0 can be seen as relative for time point 0 of the track. For new cells produced from mitosis at later time points, they are inserted into the list TP_coord_meta[0], so they can be lineage traced [FROM THEIR OWN STARTING TIME POINT, ENCODED IN 4D COORDINATE]
		####Append 3char lineage ID (for family tracing) and a gen ID... TO add generational ID, append generation 0 to everything in base origin list (before any mitotic additions) before this for loop begins. then read for the generational coordinate deeper in the loop for when +1 needs to be added to the generation.
		window.refresh()
		tracked = a_crd #TP_coord_meta[0][point2track]
		print("family_ID_clock = ", family_ID_clock)
		if tracked[4] == 0:
			tracked[4] = family_ID_clock
		else:
			pass
		family_ID_clock += 1
		print("family_ID_clock += 1", family_ID_clock)
		#print("tracked = ", a_crd)
		b = [] ###### this is the transient list that will hold the coordinates for this particular track. Once track terminates, this is appended to master track list (the list holding all the tracks)
		b.append(tracked)
		#print("pre-looped TRACK = ", tracked)
		#print("pre-looped track 3D ONLY: ", tracked[:3])
		#print("tracked[3] +1", tracked[3]+1)
		tracked_time = tracked[3]
		print("tracked_time = ", tracked_time)
		print("len tracked =  ", len(tracked))
		for a in TP_coord_meta[int(tracked_time+1)]:
		 ####to change this to look at TP_coord_meta[(time coordinate)+1], to accound for new mitotic lineages reinserted into TP_coord_meta[0], making TP_coord_meta[0] less about time point 0 in absolute terms, but time point 0 in relative terms to the lineage being traced, ergo it is an origin list. perhaps rename to some variant of origin list
			tracked_time += 1
			if tracked_time + 1 == T_len:
				print("out of time")
				break
			else:
				pass
			#print(" 'a' check, a = ", a)
			print("new tracked time = ", tracked_time)
			#print("TP_coord_meta[int(tracked[3]+1)] = ", TP_coord_meta[int(tracked[3]+1)])
			current_time_coords = TP_coord_meta[int(tracked_time+1)]

			#######################
			#scan_radius = 23-(((a[4]).count('x'))*4)  ##### SET DISTANCE TO LOOK FOR NEAREST POINTS WITHIN
			print("generation check tracked: ", tracked)
			if tracked[5]>cycles_to_plataeu: #scan radius needs to reduce with succesive cell divisions, this is dataset dependent. in this dataset we begin from 128 cell. normally we anticipate the cell diameter to stabilise past 1000 cell stage, so up until around 512 cell we reduce diameter by  approx 20percent. Past 1000  cell stage there is less need to reduce diameter so drastically
				scan_radius=scan_rad_lim #radius needn't go lower
			else:
				scan_radius=scan_rad_start-(tracked[5]*scan_rad_decay_factor) #tracked[5] contains generation information. this particular paramters will have to be revisited per dataset. depends on which generations we are looking at and the degree of cell size reduction over time. will need to inform with literature, or maybe construct class to parse cell membrane markers, if  any are availaible. this would need to be specified as a parameter.
			#scan_radius=15
			#mitosis_buffer = 15

			#print("a = ", a)

			#print("tracked[:3] = ", tracked[:3])
			#print("np.array(tracked[:3]) = ", np.array(tracked[:3]))

			#########
			print("tracked = ", tracked)
			tracked_remain = tracked[3:] #### 2019/10/17 tracked_remain remembes the tracking-irrelevant coordinates and stores them to append back onto the array when it is reappliued to master track list (so we retain generational information)
			print("tracked_remain = ", tracked[3:])
			tracked_3D_coord = SLICE(tracked)
			print("tracked 3d coord = ",tracked_3D_coord)
			current_time_coords = np.array(current_time_coords)
			#print("current_time_coords[:,:3] = ", current_time_coords[:,:3])

			numnearpts = num_pts_within_dist(tracked_3D_coord, current_time_coords[:,:3], scan_radius)
			IDnearpts = ID_pts_within_dist(tracked_3D_coord, current_time_coords[:,:3], scan_radius)
			print("IDnearpts = ", IDnearpts)
			print("LEN IDnearpts = ", len(IDnearpts))
			print("len(b)", len(b))
			if len(IDnearpts) == 1 or ((len(b) <= mitosis_buffer) and len(IDnearpts) > 0): ###### 2019-11-18 if only point is located in scan radius, proceed finding next point as usual. ALSO, newly added, the number of points needs to be above 0, otherwise the track is terminated (in else statement). If there is more than 1 item in radius but mitosis has just occurred, then proceed as if only item is nearby. This is to prevent newly divided daughters from detecting each other and registering another mitotic event.
				print("single hit in radius")
				tracked  = nearest_K_pts_v2(SLICE(tracked),np.array(IDnearpts),1)
				print("current_time_coords[0][3] =", current_time_coords[0][3])
				print("tracked = ", tracked)
				tracked = STICK(tracked, tracked_remain)
				tracked[3] = tracked_remain[0]+1
				print("tracked time + 1", tracked)
				b.append(tracked)
				if current_time_coords[0][3] + 1 == T_len:
					print("reached end of time lapse")
					break
				else:
					print("more timelapse to go")
					pass
				pass

			elif len(IDnearpts) > 1 and len(b) > mitosis_buffer:
				print(">1 points found")
				#find nearest point a
				#scan radius of a
				#if 1 other point in radius, mitosis. (record data of current point and thenear point) and carry on
				#If no, move to next nearest and try again.
				####
				####ORGANISE NEIGHBOURS IN LIST BY PROXIMITY.#### THIS IS DONE AT THE FUNCTION LEVEL IN ID_pts_within_dist using return_dist AND sort_results arguments.
				####
				list_limit = 0 #only parse closest 3 neighbours
				while list_limit<2: #only parse closest 3 neighbours
					for neighbour in IDnearpts: #only parse closest 3 neighbours
						mitosis_trigger = 0 # 2020-12-11 set up a trigger to inform the if statement which continues tracking (the false positive block). If mitosis is confirmed, then this value changes to 1, which bypases the false statement block.
						list_limit+=1
						IDnearptsN = ID_pts_within_dist(SLICE(neighbour), current_time_coords[:,:3], scan_radius/2)
						if len(IDnearptsN) == 1: #mitosis confirmation
							IDnearptsN = [neighbour, IDnearptsN[0]]
							tracked_prepend = tracked
							tracked_prepend[5] = tracked_prepend[5]+1
							daughterID = 1
							for R in IDnearptsN:
								print("checking tracked value...", tracked)
								print("R = ", R)
								#R = R.tolist()
								print("tolist R = ", R)
								R = STICK(R, tracked_remain)
								print("R STICK = ", R)
								R[3] = R[3] + 1

								#R[4] = int(str(R[4]) + str('000') + str(daughterID))
								R[4] = str(R[4]) + str('x') + str(daughterID)

								R[5] = R[5] + 1 #### +1 to gen ID
								print("final R = ", R)
								TP_coord_meta[0].append(R)
								daughterID += 1
								mitosis_trigger=1
								#print("check TP_coord_meta[0]", TP_coord_meta[0])
							#### append family ID of previous coordinate to these new tracks.
							#### add in sibling ID as a part of lineage ID
							print("end of this track")
							break

				if mitosis_trigger!=1: #false positive block
					print("continue with track, trigger was noise.")
					tracked  = nearest_K_pts_v2(SLICE(tracked),np.array(IDnearpts),1)
					print("current_time_coords[0][3] =", current_time_coords[0][3])
					print("tracked = ", tracked)
					tracked = STICK(tracked, tracked_remain)
					tracked[3] = tracked_remain[0]+1
					print("tracked time + 1", tracked)
					b.append(tracked)
				elif mitosis_trigger==1:
					print("mitosis event registered, end this track")
					break
				if current_time_coords[0][3] + 1 == T_len:
					print("reached end of time lapse")
					break
				else:
					print("more timelapse to go")
					pass
				pass

			else:
				print("no track")
				break


		track_list.append(b)
		print("len(track_list) = ", len(track_list))
		#print("track_list = ", track_list[:])
	#print("TP_coord_meta[0] final tally = ", TP_coord_meta[0])
	os.chdir(trk_output)

	with open("track_list.csv", "w", newline="") as f:
		trackwrite = csv.writer(f)

		trackwrite.writerow(["FAM_ID", "GEN", "LEN"])
		for i in track_list:
			print("csv write, 'i' = ", i)
			trackwrite.writerow([i[0][4],i[0][5],len(i)])

	restack_master = []
	for J in track_list:
		#print("J = ", J)

		#coord = narray(J)
		coord=J

		#print(coord)

		#coordA = coord.astype(np.float)
		coordA=coord

		#print(coordA[0])
		print("print coordA[0][0] = ", coordA[0][0])
		x_list = []
		y_list = []
		z_list = []
		t_list = []
		gen_list = []
		HL3 = 0
		while HL3 < len(J):
			print("len(J) = ", len(J))
			for e in coordA:
				print('e = ', e)
				xcoord = e[0]
				ycoord = e[1]
				zcoord = e[2]
				tcoord = e[3]
				gen_coord = e[5]
				print(xcoord)
				print(ycoord)
				print(zcoord)
				print(gen_coord)
				x_list.append(xcoord)
				y_list.append(ycoord)
				z_list.append(zcoord)
				t_list.append(tcoord)
				gen_list.append(gen_coord)
				print('x_list = ', x_list)
				print('y_list = ', y_list)
				print('z_list = ', z_list)
				print('gen_list = ', gen_list)
				restack_list = []
				restack_list.append(x_list) #restack_master[K][0]
				restack_list.append(y_list) #restack_master[K][1]
				restack_list.append(z_list) #restack_master[K][2]
				restack_list.append(t_list) #restack_master[K][3]
				restack_list.append(gen_list) #restack_master[K][4]
				print('#####list check: ', restack_list)
				HL3+=1
		restack_master.append(restack_list)
	#print("len track list: ", len(track_list))
	#print("###RESTACK_MASTER =", restack_master, "#########")
	window.refresh()
	fig = plt.figure() ###########################XXXXXXXXX###########################
	ax = plt.axes(projection = '3d')

	with open("tracks_.csv","w", newline="") as f:
		trackwrite = csv.writer(f)
		trackwrite.writerow(["X","Y","Z","T","GEN"])
		for K in range(0,len(restack_master)):
			print("RM[K][0] = ", restack_master[K][0])
			print("RM[K][1] = ", restack_master[K][1])
			print("xyz = ", restack_master[K][0],restack_master[K][1], restack_master[K][2])
			print("[RM[K][3][0] = ", restack_master[K][4][0])
			ax.plot(restack_master[K][0], restack_master[K][1], restack_master[K][2],c=cm.tab20c((int(restack_master[K][4][0]))))
			#trackwrite.writerow(["X","Y","Z"])
			trackwrite.writerow([restack_master[K][0], restack_master[K][1], restack_master[K][2], restack_master[K][3], restack_master[K][4]]) #contains coordinates,  time and generation for each track

	#line_ani = animation.FuncAnimation(fig, update_lines, 10, interval=50, blit=False)

	end = time.time()
	print("execution time = ", end - start,"s")
	window.refresh()
	window.read(50)
	plt.show(block=False) #2020-11-18 added block=false to synergise with GUI multi window/multi threading functionality. We do not want the matplotlib window to stall the GUI.
	#plt.show()











def fetchcolvalue(loaded_csv, row_name):
	file = open(loaded_csv, 'r')
	for row in csv.reader(file):
		if row[0] == row_name:
			return row[1]

#####
#autoload save state


loaded_csv = "saturnia_sav.csv"
try:
	if path.exists(loaded_csv):
		#try:
		input, output, max_nuc_diam, min_bright = fetchcolvalue(loaded_csv,"input"),fetchcolvalue(loaded_csv,"output"),fetchcolvalue(loaded_csv,"max_nuc_diam"),fetchcolvalue(loaded_csv,"min_bright")

		prom, peak_dist, min2Darea, max2Darea, cores = fetchcolvalue(loaded_csv,"prom"),fetchcolvalue(loaded_csv,"peak_dist"),fetchcolvalue(loaded_csv,"min2Darea"),fetchcolvalue(loaded_csv,"max2Darea"),fetchcolvalue(loaded_csv,"cores")


		conv_input, conv_output, min_vol_filt= fetchcolvalue(loaded_csv,"conv_input"),fetchcolvalue(loaded_csv,"conv_output"),fetchcolvalue(loaded_csv,"min_vol_filt")

		x_conv, y_conv, z_conv = fetchcolvalue(loaded_csv,"x_conv"),fetchcolvalue(loaded_csv,"y_conv"),fetchcolvalue(loaded_csv,"z_conv")


		trk_input, trk_output, frame_no= fetchcolvalue(loaded_csv,"trk_input"),fetchcolvalue(loaded_csv,"trk_output"),fetchcolvalue(loaded_csv,"frame_no")

		scan_rad, min_scanR, decay_rate, mit_buff = fetchcolvalue(loaded_csv,"scan_rad"),fetchcolvalue(loaded_csv,"min_scanR"),fetchcolvalue(loaded_csv,"decay_rate"),fetchcolvalue(loaded_csv,"mit_buff")
	else:
		input, output, max_nuc_diam, min_bright, prom, peak_dist, min2Darea, max2Darea, cores = "","","","","","","","",""
		conv_input, conv_output, min_vol_filt, x_conv, y_conv, z_conv = "","","","","",""
		trk_input, trk_output, frame_no, scan_rad, min_scanR, decay_rate, mit_buff = "","","","","","",""
except:
	input, output, max_nuc_diam, min_bright, prom, peak_dist, min2Darea, max2Darea, cores = "","","","","","","","",""
	conv_input, conv_output, min_vol_filt, x_conv, y_conv, z_conv = "","","","","",""
	trk_input, trk_output, frame_no, scan_rad, min_scanR, decay_rate, mit_buff = "","","","","","",""
	
#####
#TAB 1: SEGMENTATION
selection_col = [[pg.Text("Input   "), pg.In(input, size=(25,2),enable_events=True, key = "-INPUT-"), pg.FolderBrowse()],
[pg.Text("Output"), pg.In(output, size=(25,2),enable_events=True, key = "-OUTPUT-"),pg.FolderBrowse()]
]


seg_params_col = [

[ pg.Text("max nuclear diameter (pixels): "), pg.In(max_nuc_diam,size=(10,1), enable_events=True, key = "-MAXNUCDIAM-")],
[ pg.Text("Minimal brightness: "),pg.In(min_bright,size=(10,1), enable_events=True, key = "-MINBRIGHT-")],
[ pg.Text("Peak prominence: "),pg.In(prom,size=(10,1), enable_events=True, key = "-PROM-")],
[ pg.Text("clustering distance (pixels): "),pg.In(peak_dist,size=(10,1), enable_events=True, key = "-PEAKDIST-")],
[ pg.Text("Min 2D area (pixels): "),pg.In(min2Darea, size=(10,1), enable_events=True, key = "-MIN2DA-")],
[ pg.Text("Max 2D area (pixels): "),pg.In(max2Darea, size=(10,1), enable_events=True, key = "-MAX2DA-")],
[ pg.Text("Processing cores: "),pg.In(cores, size=(10,1), enable_events=True, key = "-CORES-")],

[pg.Button("Apply all", enable_events=True, key = "apply_all"),
pg.Button("Run segmentation", enable_events=True, key = "-SEGMENTATION_2D-")],

[pg.Input(visible=False, enable_events=True, key='-SAVE_SEG-'), pg.FileSaveAs("Save seg settings", enable_events=True,key="-SAVE_SEG-", file_types=(("CSV","*.csv"))),
pg.Input(visible=False, enable_events=True, key='-LOAD_SEG-'), pg.FileBrowse("Load seg settings",key="-LOAD_SEG-")
]

]

#####
#TAB 2: MATLAB MESSAGE

matlab_msg = [[pg.Text("Run MATLAB SCRIPT on segmentation outputs:\ninput the 2D segmentation output folder\nas the 3D segmentation input (if not auto-filled already)\non line 6 of the matlab code for the variable seg_input")]]
matlab_msg2 = [[pg.Text("Run MATLAB SCRIPT on segmentation outputs:\ninput the 2D segmentation output folder\nas the 3D segmentation input (if not auto-filled already)\non line 6 of the matlab code for the variable seg_input")]]

#####
#TAB 3: TRANSFORMATION
transf_input = [[pg.Text("Input   "), pg.In(conv_input, size=(25,2),enable_events=True, key = "-T_INPUT-"), pg.FolderBrowse()],
[pg.Text("Output"), pg.In(conv_output, size=(25,2),enable_events=True, key = "-T_OUTPUT-"),pg.FolderBrowse()]
]

tran_params_col = [

[ pg.Text("Min volume filter: "),pg.In(min_vol_filt, size=(10,1), enable_events=True, key = "-MINVOLFILT-")],
[ pg.Text("X conversion (micron/px): "),pg.In(x_conv, size=(10,1), enable_events=True, key = "-XCONV-")],
[ pg.Text("Y convresion (micron/px): "),pg.In(y_conv, size=(10,1), enable_events=True, key = "-YCONV-")],
[ pg.Text("Z conversion (micron/px): "),pg.In(z_conv, size=(10,1), enable_events=True, key = "-ZCONV-")],

[pg.Button("Apply all", enable_events=True, key = "apply_all_conv"),
pg.Button("Transform features", enable_events=True, key="-TRANSFORM-")],

[pg.Input(visible=False, enable_events=True, key='-SAVE_TRAN-'), pg.FileSaveAs("Save tran settings", enable_events=True,key="-SAVE_TRAN-", file_types=(('CSV','.csv'))),
pg.Input(visible=False, enable_events=True, key='-LOAD_TRAN-'), pg.FileBrowse("Load tran settings",key="-LOAD_TRAN-")
]
]

#####
#TAB 4: TRACKING
track_input = [[pg.Text("Input   "), pg.In(trk_input, size=(25,2),enable_events=True, key = "-TRK_INPUT-"), pg.FileBrowse()],
[pg.Text("Output"), pg.In(trk_output, size=(25,2),enable_events=True, key = "-TRK_OUTPUT-"),pg.FolderBrowse()]
]

trax_params_col = [

[ pg.Text("Number of frames: "),pg.In(frame_no, size=(10,1), enable_events=True, key = "-DATALEN-")],
[ pg.Text("Scan radius: "),pg.In(scan_rad, size=(10,1), enable_events=True, key = "-SCANR-")],
[ pg.Text("Min scan radius: "),pg.In(min_scanR, size=(10,1), enable_events=True, key = "-MINSCANR-")],
[ pg.Text("Decay factor: "),pg.In(decay_rate, size=(10,1), enable_events=True, key = "-SCANDECAY-")],
[ pg.Text("Mitosis buffer: "),pg.In(mit_buff, size=(10,1), enable_events=True, key = "-MITBUFF-")],
[pg.Button("Apply all", enable_events=True, key = "apply_all_track"),
pg.Button("Track", enable_events=True, key="-TRACK_LIN-")],

[pg.Input(visible=False, enable_events=True, key='-SAVE_TRACK-'), pg.FileSaveAs("Save track settings", enable_events=True,key="-SAVE_TRACK-", file_types=(('CSV','*.csv'))),
pg.Input(visible=False, enable_events=True, key='-LOAD_TRACK-'), pg.FileBrowse("Load track settings", key="-LOAD_TRACK-")
]

]
#####
#CONSOLE OUTPUT
image_view_col = [[pg.Text("Console Output"), pg.Input(visible=False, enable_events=True, key='-SAVE_ALL-'), pg.FileSaveAs("Save all settings", enable_events=True,key="-SAVE_ALL-", file_types=(('CSV','*.csv'))),
pg.Input(visible=False, enable_events=True, key='-LOAD_ALL-'), pg.FileBrowse("Load all settings", key="-LOAD_ALL-"),
pg.Checkbox("save all settings on exit",enable_events=True, key="exit_save")],
[pg.Output(size = (50,10), key="-CONSOLE1-")],
[pg.Button("Clear all fields", enable_events=True, key="-CLEAR_ALL-"),pg.Button("Clear console")]]



layout_tab1=[[pg.Column(selection_col), pg.VSeperator(),pg.Column(seg_params_col)]]
layout_tab2=[[pg.Column(matlab_msg), pg.VSeperator(), pg.Column(matlab_msg2)]]
layout_tab3=[[pg.Column(transf_input),pg.VSeperator(),pg.Column(tran_params_col)]]
layout_tab4=[[pg.Column(track_input),pg.VSeperator(),pg.Column(trax_params_col)]]
#layout_console_1 = [[pg.Column(image_view_col)]]

layout = [[pg.TabGroup([[pg.Tab("2D Segmentation", layout_tab1),pg.Tab("3D Segmentation", layout_tab2),pg.Tab("Transformation", layout_tab3),pg.Tab("Tracking", layout_tab4)]])
,pg.Column(image_view_col)]]


window=pg.Window("Saturnia",layout)

def value_check(valname, value):
	try:
		value = float(value)
		if value%1==0:
			value=int(value)
	except:
		value = value
	if isinstance(value, int) or isinstance(value,float):
		print("{} : {} selected".format(valname, value))
		return value
	else:
		print("{} PARAM NOT SELECTED".format(valname))




def mod_run_chk(inP, outP, var_list, run_msg, all_systems):
	#segmentation_vars = var_list
	all_systems=""
	segment_run_status=""
	for segment_var in var_list:
		try:
			segment_var = float(segment_var)
			if segment_var%1==0:
				segment_var=int(segment_var)
		except:
			segment_var = segment_var
		if segment_var=="" or isinstance(segment_var,str):
			print("type problem ", segment_var)
			print(type(segment_var))
			segment_run_status = "NOT_READY"
			#print(segment_run_status)
	'''
	if inP=="" or outP=="":
		segment_run_status = "NOT_READY"
	'''
	if path.exists(inP) and path.exists(outP):
		pass
	else:
		segment_run_status = "NOT_READY"
		print(inP)
		print(outP)
		print("IO path error")
	if segment_run_status=="NOT_READY":
		print("Missing input parameter(s)")
	else:
		print(run_msg)
		#a_s = "GO"
		all_systems = "GO"
		#print("post run msg")
		#print("function all_systems = ", all_systems)
		return all_systems
			#mod_run_chk(input,output,[max_nuc_diam, min_bright, prom, peak_dist, min2Darea, max2Darea],all_systems)

while True:
	event, values = window.read(20)
	#event == "Exit"
	if event == pg.WIN_CLOSED:

		all_file_basename = "saturnia_sav.csv"
		all_val_list = [input,output,max_nuc_diam,min_bright,prom,peak_dist,min2Darea,max2Darea,cores,
		conv_input,conv_output,min_vol_filt,x_conv,y_conv,z_conv,
		trk_input,trk_output,frame_no,scan_rad,min_scanR,decay_rate,mit_buff
		]

		all_valName_list = ["input","output","max_nuc_diam","min_bright","prom","peak_dist","min2Darea","max2Darea","cores",
		"conv_input","conv_output","min_vol_filt","x_conv","y_conv","z_conv",
		"trk_input","trk_output","frame_no","scan_rad","min_scanR","decay_rate","mit_buff"
		]

		with open(all_file_basename, "w") as csv_file:
			writer = csv.writer(csv_file,delimiter=',')
			for idx, row in enumerate(all_val_list):
				print(all_valName_list[idx],row)
				writer.writerow([all_valName_list[idx],row])
				#csv_file.write(seg_valName_list[idx],row)
				#writer.writerow(row)
		csv_file.close()
		#os.chdir(gui_curDir)

		break

	##SEGMENTATION##
	if event == "-INPUT-":
		input=values["-INPUT-"]+"/"
		print("segment input folder= ", input)
	elif event == "-OUTPUT-":
		output=values["-OUTPUT-"]+"/"
		print("segment output folder = ", output)
	if event == "apply_all":
		max_nuc_diam = values["-MAXNUCDIAM-"]
		value_check("Max nuclear diameter", max_nuc_diam)
		min_bright = values["-MINBRIGHT-"]
		value_check("minimum brightness", min_bright)
		prom = values["-PROM-"]
		value_check("peak prominence", prom)
		peak_dist = values["-PEAKDIST-"]
		value_check("distance between peaks", peak_dist)
		min2Darea = values["-MIN2DA-"]
		value_check("Minimum 2D area", min2Darea)
		max2Darea = values["-MAX2DA-"]
		value_check("Maximum 2D area", max2Darea)
		cores = values["-CORES-"]
		value_check("cores", cores)
		#insert segmentation code reference here
		print("___________\n___________\n")

	if event=="-SAVE_SEG-":
		if not values["-SAVE_SEG-"]:
			print("save seg cancelled")
		elif values["-SAVE_SEG-"]:
			seg_file_basename = values["-SAVE_SEG-"]
			seg_val_list = [input,output,max_nuc_diam,min_bright,prom,peak_dist,min2Darea,max2Darea,cores]
			seg_valName_list = ["input","output","max_nuc_diam","min_bright","prom","peak_dist","min2Darea","max2Darea","cores"]

			with open(seg_file_basename, "w") as csv_file:
				writer = csv.writer(csv_file,delimiter=',')
				for idx, row in enumerate(seg_val_list):
					print(seg_valName_list[idx],row)
					writer.writerow([seg_valName_list[idx],row])
					#csv_file.write(seg_valName_list[idx],row)
					#writer.writerow(row)
			csv_file.close()
			#os.chdir(gui_curDir)
		else:
			print("oh noe")

	if event == "-LOAD_SEG-" or event =="load_all":
		loaded_csv = values["-LOAD_SEG-"]
		input=fetchcolvalue(loaded_csv,"input")
		window["-INPUT-"].update(input)
		output=fetchcolvalue(loaded_csv,"output")
		window["-OUTPUT-"].update(output)
		max_nuc_diam=fetchcolvalue(loaded_csv, "max_nuc_diam")
		window["-MAXNUCDIAM-"].update(max_nuc_diam)
		min_bright=fetchcolvalue(loaded_csv, "min_bright")
		window["-MINBRIGHT-"].update(min_bright)
		prom=fetchcolvalue(loaded_csv, "prom")
		window["-PROM-"].update(prom)
		peak_dist=fetchcolvalue(loaded_csv, "peak_dist")
		window["-PEAKDIST-"].update(peak_dist)
		min2Darea=fetchcolvalue(loaded_csv, "min2Darea")
		window["-MIN2DA-"].update(min2Darea)
		max2Darea=fetchcolvalue(loaded_csv, "max2Darea")
		window["-MAX2DA-"].update(max2Darea)
		cores=fetchcolvalue(loaded_csv, "cores")
		window["-CORES-"].update(cores)

	if event == "-SEGMENTATION_2D-":
		all__systems="" #need to declare this for the mod_run_chk to work, as it is the indicator variable that the imported module can run
		#print("sanity check")
		all__systems = mod_run_chk(input,output,[max_nuc_diam, min_bright, prom, peak_dist, min2Darea, max2Darea,cores],"running!",all__systems)
		print("--all_systems= ", all__systems)
		if all__systems=="GO":
			print("running!!!!!")
			try:
				gui_dir = os.getcwd()
				input_dir=input
				output_dir=output
				print("segmentation_")
				curdirA=input_dir
				os.chdir(curdirA)
				print("curdirA = ", curdirA)
				print("list = ", os.listdir(curdirA))
				for dir in os.listdir(curdirA):
					try:
						print(os.path.abspath(dir))
						dirb = os.path.abspath(dir)
						#print(dirb)
						print(dir)
						if os.path.isdir(dir) is True:
							print("XXY")
							output_dir = output+"SEG"+dir+"/" #using output rather than outputdir as output represents initally stated path for output, and output_dir gets changed, leading to nesting of directories in the first directory scanned.
							if path.exists(output_dir):
								output_copy = 0
								print("exist_1")
								while path.exists(output_dir):
									output_dir=output_dir[:-1]+string(output_copy)+"/"
									print("exist_2")
									output_copy+=1
							print("XXY")
							print(output_dir)
							os.mkdir(output_dir)
							print("XXY")
							curdir=os.chdir(dir)
							print("XXY")
							listdirec=os.listdir(curdir)
							print("XXYZ")
							#with Pool(cores) as p:


								#p.map(sliceBySlice(a,output_dir,max_nuc_diam,min_bright,peak_dist,prom,min2Darea,max2Darea,cores),listdirec)
								#os.chdir(curdirA)

							'''
							jerbs = []
							#que = Queue()
							def initiate_process(procesX):
								procesX.start()
							'''
							input_list = [f'{str(input_dir+a)},{output_dir},{max_nuc_diam},{min_bright},{peak_dist},{prom},{min2Darea},{max2Darea},{cores}'.split(",") for a in listdirec]

							#print("INPUT_LIST=",input_list)

							#print("input list = ", input_list)
							try:
									#with Pool(processes=int(cores),maxtasksperchild=500) as p:

								with closing(Pool(int(cores))) as p: #close processes and recover memory once each process is done
									p.starmap(sliceBySlice,input_list)
									window.refresh()
									#p.close()
									#p.join()
									#window.read(10)
									#window.refresh()
									os.chdir(curdirA)


								#pool.close()

								'''
								p =mp.Pool(2)
								#for a in range(2):
								#print(a)
								a = input_list
								window.refresh()
								print(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7])
								#window.refresh()
								p.map_async(sliceBySlice, args=(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8]))

								#p.close()
								#p.join()
								window.refresh()

								os.chdir(curdirA)
								'''
							except Exception as exception:
								# Get current system exception
								ex_type, ex_value, ex_traceback = sys.exc_info()

								# Extract unformatter stack traces as tuples
								trace_back = traceback.extract_tb(ex_traceback)

								# Format stacktrace
								stack_trace = list()

								for trace in trace_back:
									stack_trace.append("File : %s , Line : %d, Func.Name : %s, Message : %s" % (trace[0], trace[1], trace[2], trace[3]))

								print("Exception type : %s " % ex_type.__name__)
								print("Exception message : %s" %ex_value)
								print("Stack trace : %s" %stack_trace)

								print("error?")

					except:
						# Get current system exception
						ex_type, ex_value, ex_traceback = sys.exc_info()

						# Extract unformatter stack traces as tuples
						trace_back = traceback.extract_tb(ex_traceback)

						# Format stacktrace
						stack_trace = list()

						for trace in trace_back:
							stack_trace.append("File : %s , Line : %d, Func.Name : %s, Message : %s" % (trace[0], trace[1], trace[2], trace[3]))

						print("Exception type : %s " % ex_type.__name__)
						print("Exception message : %s" %ex_value)
						print("Stack trace : %s" %stack_trace)

						print("error?")
						#print("dir error")
						#func_thread_seg = threading.Thread(target=sliceBySlice, args=(input,output,max_nuc_diam,min_bright,peak_dist,prom,min2Darea,max2Darea,cores),daemon=True)
						#func_thread_seg.start()
						#^^^ pysimpleguiQt needs initiate parallel threads to run new scripts, as it must remain the primary thread.
						#func_thread_seg.join()
					print("post seg")
						#window.refresh()

					#sp=subprocess.Popen(segment_2D_stack(input,output,max_nuc_diam,min_bright,peak_dist,prom,min2Darea,max2Darea,cores),shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=subprocess.PIPE)
					#segment_2D_stack(input,output,max_nuc_diam,min_bright,peak_dist,prom,min2Darea,max2Darea,cores)
			except:
				# Get current system exception
				ex_type, ex_value, ex_traceback = sys.exc_info()

				# Extract unformatter stack traces as tuples
				trace_back = traceback.extract_tb(ex_traceback)

				# Format stacktrace
				stack_trace = list()

				for trace in trace_back:
					stack_trace.append("File : %s , Line : %d, Func.Name : %s, Message : %s" % (trace[0], trace[1], trace[2], trace[3]))

				print("Exception type : %s " % ex_type.__name__)
				print("Exception message : %s" %ex_value)
				print("Stack trace : %s" %stack_trace)
				print("class error?")
			os.chdir(gui_dir)



	##TRANSFORMATION##
	if event == "-T_INPUT-":
		conv_input=values["-T_INPUT-"]+"/"
		print("transform input folder= ", conv_input)
		conv_output = conv_input
		window["-T_OUTPUT-"].update(conv_output)
	elif event == "-T_OUTPUT-":
		if values["-T_OUTPUT-"].endswith("/"):
			conv_output=values["-T_OUTPUT-"]
		else:
			conv_output=values["-T_OUTPUT-"]+"/"
		print("transform output folder = ", conv_output)
	if event == "apply_all_conv":
		min_vol_filt = values["-MINVOLFILT-"]
		value_check("Minimum volume", min_vol_filt)
		x_conv = values["-XCONV-"]
		value_check("x axis conversion", x_conv)
		y_conv = values["-YCONV-"]
		value_check("y axis conversion", y_conv)
		z_conv = values["-ZCONV-"]
		value_check("z axis conversion", z_conv)
		print("___________\n___________\n")

	if event=="-SAVE_TRAN-":
		if not values["-SAVE_TRAN-"]:
			print("save transformation params cancelled")
		elif values["-SAVE_TRAN-"]:
			tran_file_basename = values["-SAVE_TRAN-"]
			tran_val_list = [conv_input,conv_output,min_vol_filt,x_conv,y_conv,z_conv]
			tran_valName_list = ["conv_input","conv_output","min_vol_filt","x_conv","y_conv","z_conv"]

			with open(tran_file_basename, "w") as csv_file:
				writer = csv.writer(csv_file,delimiter=',')
				for idx, row in enumerate(tran_val_list):
					print(tran_valName_list[idx],row)
					writer.writerow([tran_valName_list[idx],row])
			csv_file.close()
		else:
			print("oh noe")

	if event == "-LOAD_TRAN-" or event =="load_all":
		loaded_csv = values["-LOAD_TRAN-"]
		conv_input=fetchcolvalue(loaded_csv,"conv_input")
		window["-T_INPUT-"].update(conv_input)
		conv_output=fetchcolvalue(loaded_csv,"conv_output")
		window["-T_OUTPUT-"].update(conv_output)
		min_vol_filt=fetchcolvalue(loaded_csv, "min_vol_filt")
		window["-MINVOLFILT-"].update(min_vol_filt)
		x_conv=fetchcolvalue(loaded_csv, "x_conv")
		window["-XCONV-"].update(x_conv)
		y_conv=fetchcolvalue(loaded_csv, "y_conv")
		window["-YCONV-"].update(y_conv)
		z_conv=fetchcolvalue(loaded_csv, "z_conv")
		window["-ZCONV-"].update(z_conv)

	if event == "-TRANSFORM-":
		all_systems="" #need to declare this for the mod_run_chk to work, as it is the indicator variable that the imported module can run
		all_systems=mod_run_chk(conv_input,conv_output,[min_vol_filt,x_conv,y_conv,z_conv],"running!",all_systems)
		if all_systems=="GO":
			try:
				outpath=""
				#feature_transform(conv_input,conv_output,min_vol_filt,x_conv,y_conv,z_conv,outpath)

				print("transformation_")
				func_thread_trnf = threading.Thread(target=feature_transform, args=(conv_input,conv_output,min_vol_filt,x_conv,y_conv,z_conv,outpath),daemon=True)
				func_thread_trnf.start()
				#^^^ pysimpleguiQt needs initiate parallel threads to run new scripts, as it must remain the primary thread.
				func_thread_trnf.join()
				print("post trnsf")

				#sp=subprocess.Popen(feature_transform(conv_input,conv_output,min_vol_filt,x_conv,y_conv,z_conv,outpath),shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=subprocess.PIPE)
				window["-TRK_INPUT-"].update(outpath)
				trk_input = outpath
			except:
					# Get current system exception
					ex_type, ex_value, ex_traceback = sys.exc_info()

					# Extract unformatter stack traces as tuples
					trace_back = traceback.extract_tb(ex_traceback)

					# Format stacktrace
					stack_trace = list()

					for trace in trace_back:
					    stack_trace.append("File : %s , Line : %d, Func.Name : %s, Message : %s" % (trace[0], trace[1], trace[2], trace[3]))

					print("Exception type : %s " % ex_type.__name__)
					print("Exception message : %s" %ex_value)
					print("Stack trace : %s" %stack_trace)

					print("error?")



	##TRACKING##
	if event == "-TRK_INPUT-":
		trk_input=values["-TRK_INPUT-"]
		print(" track input= ", trk_input)
	elif event == "-TRK_OUTPUT-":
		trk_output=values["-TRK_OUTPUT-"]+"/"
		print(" track output folder = ", trk_output)
	if event == "apply_all_track":
		frame_no = values["-DATALEN-"]
		value_check("frames", frame_no)
		scan_rad = values["-SCANR-"]
		value_check("scan radius", scan_rad)
		min_scanR = values["-MINSCANR-"]
		value_check("Minimum scan radius", min_scanR)
		decay_rate = values["-SCANDECAY-"]
		value_check("rate of decay", decay_rate)
		mit_buff = values["-MITBUFF-"]
		value_check("Mitosis time buffer", mit_buff)
		print("___________\n___________\n")
	if event=="-SAVE_TRACK-":
		if not values["-SAVE_TRACK-"]:
			print("save track params cancelled")
		elif values["-SAVE_TRACK-"]:
			track_file_basename = values["-SAVE_TRACK-"]
			track_val_list = [trk_input,trk_output,frame_no,scan_rad,min_scanR,decay_rate,mit_buff]
			track_valName_list = ["trk_input","trk_output","frame_no","scan_rad","min_scanR","decay_rate","mit_buff"]

			with open(track_file_basename, "w") as csv_file:
				writer = csv.writer(csv_file,delimiter=',')
				for idx, row in enumerate(track_val_list):
					print(track_valName_list[idx],row)
					writer.writerow([track_valName_list[idx],row])
			csv_file.close()
		else:
			print("oh noe")

	if event == "-LOAD_TRACK-" or event =="load_all":
		loaded_csv = values["-LOAD_TRACK-"]
		trk_input=fetchcolvalue(loaded_csv,"trk_input")
		window["-TRK_INPUT-"].update(trk_input)
		trk_output=fetchcolvalue(loaded_csv,"trk_output")
		window["-TRK_OUTPUT-"].update(trk_output)
		frame_no=fetchcolvalue(loaded_csv, "frame_no")
		window["-DATALEN-"].update(frame_no)
		scan_rad=fetchcolvalue(loaded_csv, "scan_rad")
		window["-SCANR-"].update(scan_rad)
		min_scanR=fetchcolvalue(loaded_csv, "min_scanR")
		window["-MINSCANR-"].update(min_scanR)
		decay_rate=fetchcolvalue(loaded_csv, "decay_rate")
		window["-SCANDECAY-"].update(decay_rate)
		mit_buff=fetchcolvalue(loaded_csv, "mit_buff")
		window["-MITBUFF-"].update(mit_buff)

	if event == "-TRACK_LIN-":
		gui_dir = os.getcwd()
		q=queue.Queue()
		all_systems="" #need to declare this for the mod_run_chk to work, as it is the indicator variable that the imported module can run
		all_systems=mod_run_chk(trk_input,trk_output,[frame_no, scan_rad, min_scanR, decay_rate, mit_buff],"running!",all_systems)
		if all_systems=="GO":
			print("tracking initiated")
			try:
				#window.read(20)
				window.refresh()
				print("meep")

				print("track_")
				#func_thread_trk = threading.Thread(target=track_nuc,args=(trk_input,trk_output,frame_no, scan_rad, min_scanR, decay_rate, mit_buff),daemon=True)
				
				#func_thread_trk = threading_norm.start_new_thread(track_nuc,(trk_input,trk_output,frame_no, scan_rad, min_scanR, decay_rate, mit_buff))
				'''
				pool = mpg.Pool(1)
				pool.apply(track_nuc,(trk_input,trk_output,frame_no, scan_rad, min_scanR, decay_rate, mit_buff))
				'''

				with closing(mpg.Pool(1)) as pl: #close processes and recover memory once each process is done
					window.read(20)
					pl.apply(track_nuc,(trk_input,trk_output,frame_no, scan_rad, min_scanR, decay_rate, mit_buff))
					pl.close()
					window.refresh()

				#pool.join()
				#^^^ pysimpleguiQt needs initiate parallel threads to run new scripts, as it must remain the primary thread.
				#func_thread_trk.join()
				print("post trk")
				os.chdir(gui_dir)


			except:
					# Get current system exception
					ex_type, ex_value, ex_traceback = sys.exc_info()

					# Extract unformatter stack traces as tuples
					trace_back = traceback.extract_tb(ex_traceback)

					# Format stacktrace
					stack_trace = list()

					for trace in trace_back:
					    stack_trace.append("File : %s , Line : %d, Func.Name : %s, Message : %s" % (trace[0], trace[1], trace[2], trace[3]))

					print("Exception type : %s " % ex_type.__name__)
					print("Exception message : %s" %ex_value)
					print("Stack trace : %s" %stack_trace)

					print("error?")


	if event == "-CLEAR_ALL-":
		input, output, max_nuc_diam, min_bright, prom, peak_dist, min2Darea, max2Darea,cores = "","","","","","","","",""
		conv_input, conv_output, min_vol_filt, x_conv, y_conv, z_conv = "","","","","",""
		trk_input, trk_output, frame_no, scan_rad, min_scanR, decay_rate, mit_buff = "","","","","","",""
		window["-INPUT-"].update(input)
		window["-OUTPUT-"].update(output)
		window["-MAXNUCDIAM-"].update(max_nuc_diam)
		window["-MINBRIGHT-"].update(min_bright)
		window["-PROM-"].update(prom)
		window["-PEAKDIST-"].update(peak_dist)
		window["-MIN2DA-"].update(min2Darea)
		window["-MAX2DA-"].update(max2Darea)
		window["-CORES-"].update(cores)
		window["-T_INPUT-"].update(conv_input)
		window["-T_OUTPUT-"].update(conv_output)
		window["-MINVOLFILT-"].update(min_vol_filt)
		window["-XCONV-"].update(x_conv)
		window["-YCONV-"].update(y_conv)
		window["-ZCONV-"].update(z_conv)		
		window["-TRK_INPUT-"].update(trk_input)
		window["-TRK_OUTPUT-"].update(trk_output)
		window["-DATALEN-"].update(frame_no)
		window["-SCANR-"].update(scan_rad)
		window["-MINSCANR-"].update(min_scanR)
		window["-SCANDECAY-"].update(decay_rate)
		window["-MITBUFF-"].update(mit_buff)



	if event == "Clear console":
		window["-CONSOLE1-"].update("")


window.close()
