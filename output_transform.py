import os
from os import path
import csv
import sys
#sysargv[1] is the volume threshold. volume below this is excluded in appendcsv output. sysargv[2], sysargv[3], sysargv[4] are x y z pixel measurements in microns
def feature_transform(input,output,vol_thr,xconv,yconv,zconv,outpath):
	inputarg1=float(vol_thr)
	inputarg2=float(xconv)
	inputarg3=float(yconv)
	inputarg4=float(zconv) #convert input arguments to float for mathematical manipulation
	#inputarg5=float(sys.argv[5]) #z length lower bound
	#inputarg6=float(sys.argv[6]) #z length upper bound

	def ArgX(rowitem,inputarg):
		rowitem=float(rowitem)
		conversion=str(rowitem*inputarg) #convert rowitem to float and multiply by input argument (also assumed to float, per the creation of the inputarg variables)
		return conversion

	cwd = input #get current directory
	print(cwd) #confirm current directory
	listdir=os.listdir(cwd) #save list of files in current directory, to loop through to find relevant .csv files

	for i in listdir:
		print(i)
		i_in = input+i
		i_out = output+i
		if "_t" in i and ".csv" in i: #if file name contains '_t'. this suggests its from an imageJ time point output. The .csv ensures said files are csv's as well.
			with open(i_in,"r") as file_in: #open current file in listdir to read
				with open(i_out[:-4]+"refX.csv","w") as file_out: #set file out. i[:-4] selects file name before file extension. concat with refX.csv to give unique identifier to all newly modified/created csv's. the refX is a handle we can search for when parsing modifed csv's for appending later on
					colID = i[-7:-4] #take out the time point from the file name. ensure imageJ exports files in XXX_t00X.tif format. MATLAB will create csv's with the right name type for this code to search.
					writer = csv.writer(file_out)
					rowID=1
					for row in csv.reader(file_in): #iterate through rows in csv
						if rowID==1: #if it is the first line, create the file headers. We are inserting the time point (from file name) as first column.
							writer.writerow(["T"]+["V"]+["X"]+["Y"]+["Z"]+["ulx"]+["uly"]+["ulz"]+["xl"]+["yl"]+["zl"]) #RENAME HEADERS. T IS THE ID (THE FILE TIMEPOINT), X Y Z ASSUMES THIS IS THE ORDER OF COORDS OUTPUTTED BY , ulx, uly, ul, xl, yl, zl are the coords for the bounding box. ul(xyz) are the upper left corner  and (xyz)l are the lengths away from the ul in their respective direction. main thing to look at is are dimensions ok> I am filtering features with a zl below 3,  as these may be noise or just not enough to track
						elif int(row[0])<inputarg1:
							pass
						else: #all other rows to be rewritten with only column ID + xyz centroids.
							writer.writerow([colID]+row[0:10])
						rowID+=1
		else:
			continue #not interested in other files

	if path.exists(output):
		out_file = output+"/appendcsv.csv"
		if out_file.endswith("//appendcsv.csv"):
			out_file = output+"appendcsv.csv"
	else:
		out_file = input

	with open(out_file,"w") as appended_output: #here we append all the modified new csv's (now with only columns, ID + xyz centroids. append in one big file, remove headers from all but first file. create our output file here
		listdir=os.listdir(cwd)
		appendtofile = csv.writer(appended_output)
		fileID=1 #important to set fileID so we know if first file or not. we must retain headers from first file while removing the headers from all appended tables.
		#print("listdir=",listdir)
		for i in listdir:
			i_out = output+i
			if "refX.csv" in i and "_t" in i: #using the refX.csv handle we set. the _t is just further locking the file choice to what we want though perhaps not strictly necessary
				with open(i_out,"r") as refX:
					readrefX = csv.reader(refX)
					hrowID=1 #so we know which row we are looking at within the csv we are iterating through.
					for row in readrefX:
						print(ArgX(0.5,inputarg2))
						if fileID==1 and hrowID==1: #if first file in list (earliest time point, we hope), then we want to write everything to the new csv. header and all.
							appendtofile.writerow(row)
						elif hrowID==1: #if fileID isn't 1, then we must be appending tables below another table, ergo we must check if we are on the first row, so we can not commit this to the table. we pass on this row to iterate round the loop to +1 to our hrowID, and start appending the data rows
							pass
						else:
							row[2]=ArgX(row[2],inputarg2);row[3]=ArgX(row[3],inputarg3);row[4]=ArgX(row[4],inputarg4) #change column values first, store new row in memory to pass into writerow
							appendtofile.writerow(row[0:]) #all other possibilities must be non header rows, we want these. add to table pls.
							print("file is", i)
						hrowID+=1
			else:
				print("file is:", i)
				continue
			fileID+=1 #iterate on file number.
		outpath = out_file
		return outpath
