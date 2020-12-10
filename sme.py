import PySimpleGUIQt as pg
from os import path
import os
import csv
#import get_variable_name
import numpy as np
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
from nucleus_track_3d_2 import *
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
	event, values = window.read()
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
				window.read(20)
				window.refresh()
				print("meep")

				print("track_")
				#func_thread_trk = threading.Thread(target=track_nuc,args=(trk_input,trk_output,frame_no, scan_rad, min_scanR, decay_rate, mit_buff),daemon=True)
				
				#func_thread_trk = threading_norm.start_new_thread(track_nuc,(trk_input,trk_output,frame_no, scan_rad, min_scanR, decay_rate, mit_buff))
				
				pool = mpg.Pool(1)
				pool.apply(track_nuc,(trk_input,trk_output,frame_no, scan_rad, min_scanR, decay_rate, mit_buff))
				

				
				
				print("meep2")
				window.read(20)
				window.Refresh()
				#func_thread_trk.start()
				print("meep3")
				window.Read(20)
				window.Refresh()
				print("meep4")
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
