#script to sort raw unsorted zeiss z1 lightsheet output file series (.czi). The problem is that the time points are saved to individual .czi files and are not referencable by a master file. This script puts brackets around the suffix numbers denoting their position in the time series. If such a file does not exist, you can generate a master file which will point to all these now sorted files (which have been given brackets). The masterfile simply is a .czi file which carries the same name as all the time series files (which you just sorted), but it lacks a suffix. This can be created by copying a random time point file and removing the number suffix. Opening this file in zen or ImageJ will open all of the time points in the appropriate series.

#to get this script to work, change the value of 'presuffixLen' to how many characters there are before the numerical suffix which denotes which time point a particular .czi file is referring to.

import os
from os import path

presuffixLen = 20

currentDir = os.getcwd()
directoryFileList= [files for files in os.listdir(currentDir) if ".czi" in files]
print(directoryFileList)

print(directoryFileList[0][-5])

for files in directoryFileList:
	if files[-5].isdigit:
		filesIns = files[:-4]+")"+files[-4:]
		filesIns = filesIns[:int(presuffixLen)]+"("+filesIns[int(presuffixLen):]
		print(files)
		print(filesIns)
		os.rename(files,filesIns)
