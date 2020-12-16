#IMPORTANT BEFORE YOU RUN ANY CODE. CHANGE VALUES FOR data AND scan_radius
import time
#import PySimpleGUIQt
start = time.time()
import matplotlib
import matplotlib.pyplot as plt ##!! dependency: ensure tkinter is installed, for matplotlib gui (apt-get install python3-tk)
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TK
import matplotlib.pylab as pylab
import numpy
from numpy import ndarray as npd
import csv
import pandas as pd
import scipy.spatial as spatial
from scipy.spatial.distance import cdist
from sklearn import neighbors as skn
import sys
import os


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
	####USER INPUTS:
	#track_input

	data_len = int(frame_no) #length of dataset in frames/timepoints
	scan_rad_lim = value_check_trkparam(min_scanR) #the lowest scan radius to decay to.
	scan_rad_start = value_check_trkparam(scan_rad) #radius to scan within, to begin decay from
	scan_rad_decay_factor = value_check_trkparam(decay_rate) #multiplied by the cell generation and subtracted from scan_rad_start
	mitosis_buffer = int(mit_buff)
	####END USER INPUTS



	min_scan_diam = value_check(min_scanR)
	cycles_to_plataeu = (scan_rad_start - min_scan_diam)/scan_rad_decay_factor #number of cell cycles before cell diameter doesn't reduce dramatically
	 #number of frames after mitosis that another mitosis event should not be triggered
	data_len = data_len + 1
	def narray(lst):
		lst = numpy.array(lst, dtype=numpy.float)
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
		arroo = numpy.array(arroo)
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
			if isinstance(arroo,numpy.ndarray):
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
	elif isinstance(TP_coord_meta, numpy.ndarray):
		print("TP_coord_meta is an array")
	else:
		print("not sure what TP_coord_meta is")
	'''
	family_ID_clock = 0
	for a_crd in TP_coord_meta[0]:          ####TP_coord_meta[0] is the list of coordinates in time point 0 but time point 0 can be seen as relative for time point 0 of the track. For new cells produced from mitosis at later time points, they are inserted into the list TP_coord_meta[0], so they can be lineage traced [FROM THEIR OWN STARTING TIME POINT, ENCODED IN 4D COORDINATE]
		####Append 3char lineage ID (for family tracing) and a gen ID... TO add generational ID, append generation 0 to everything in base origin list (before any mitotic additions) before this for loop begins. then read for the generational coordinate deeper in the loop for when +1 needs to be added to the generation.
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
			#print("numpy.array(tracked[:3]) = ", numpy.array(tracked[:3]))

			#########
			print("tracked = ", tracked)
			tracked_remain = tracked[3:] #### 2019/10/17 tracked_remain remembes the tracking-irrelevant coordinates and stores them to append back onto the array when it is reappliued to master track list (so we retain generational information)
			print("tracked_remain = ", tracked[3:])
			tracked_3D_coord = SLICE(tracked)
			print("tracked 3d coord = ",tracked_3D_coord)
			current_time_coords = numpy.array(current_time_coords)
			#print("current_time_coords[:,:3] = ", current_time_coords[:,:3])

			numnearpts = num_pts_within_dist(tracked_3D_coord, current_time_coords[:,:3], scan_radius)
			IDnearpts = ID_pts_within_dist(tracked_3D_coord, current_time_coords[:,:3], scan_radius)
			print("IDnearpts = ", IDnearpts)
			print("LEN IDnearpts = ", len(IDnearpts))
			print("len(b)", len(b))
			if len(IDnearpts) == 1 or ((len(b) <= mitosis_buffer) and len(IDnearpts) > 0): ###### 2019-11-18 if only point is located in scan radius, proceed finding next point as usual. ALSO, newly added, the number of points needs to be above 0, otherwise the track is terminated (in else statement). If there is more than 1 item in radius but mitosis has just occurred, then proceed as if only item is nearby. This is to prevent newly divided daughters from detecting each other and registering another mitotic event.
				print("single hit in radius")
				tracked  = nearest_K_pts_v2(SLICE(tracked),numpy.array(IDnearpts),1)
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
					tracked  = nearest_K_pts_v2(SLICE(tracked),numpy.array(IDnearpts),1)
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

		#coordA = coord.astype(numpy.float)
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
	#plt.show(block=False) #2020-11-18 added block=false to synergise with GUI multi window/multi threading functionality. We do not want the matplotlib window to stall the GUI.
	plt.show()
