import os

currDir = os.getcwd()

dirList = os.listdir(currDir)

#print(dirList)
dirList.remove('imageFolderSort.py')
prefixes = []
for file in dirList:
	if file.endswith(".tif"):
		prefixes.append(file)
		
prefixes = [w[:-9] for w in prefixes]
prefix_set = set(prefixes)
print(prefix_set)

for a in prefix_set:
	if a[:3] == '16c':
		print(a)
		os.mkdir(a)
