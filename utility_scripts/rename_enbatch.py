#for appropriately renaming misnameed lightsheet output czi files, which have been numbered incorrectly. orders them by time created and numbers them accordingly. Then a master file can be created by copying a random czi from this batch and removing its numbering. 

import os
import glob

files_=glob.glob("*.czi")
files_.sort(key=os.path.getmtime)
#print(files_)

timer = 1

for a in files_:
	#print(a)
	a_rename = a[:17]+"("+str(timer)+")"+".czi"
	print(a_rename)
	os.rename(a,a_rename)
	timer+=1
