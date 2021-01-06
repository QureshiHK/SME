import os
import shutil

curDir = os.getcwd()
print(curDir)
listDir = os.listdir(curDir)

filesI = [f for f in listDir if os.path.isfile(f)]

print(filesI)


for a in filesI:
	print(a)
	print(a[:-4])
	print(a[:len(a)-9])
	if a[-4:] == '.tif':
		shutil.move(a,str(a[:len(a)-9])+'/')
