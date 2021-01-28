#program to read lineages obtained from tracking program. requires the track_list.csv output. compare lineages and potentially draw lineage trees

import pandas as pd
import numpy as np
from numpy.polynomial.polynomial import polyfit
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt

lowlimX=65 #exclude nuclei with cell cycles which fall below this in X axis
highlimX=100 #exclude nuclei with cell cycles which exceed this in X axis
lowlimY=65 #exclude nuclei with cell cycles which fall below this in Y axis
highlimY=115 #exclude nuclei with cell cycles which exceed this in Y axis
dataset="wtPii" #sets dataset title for chart title
gen_look_1 = int(1) #generation plotted in x axis. Gen 0 is not plotted as that is considered 'incomplete' first generation.
gen_look_2 = int(2) #generation plotted in y axis
#gen_look_3 = int(3)
frameR=10 #frame rate is time in seconds between frames. this determines the units to use on axes


conversionFrameR=frameR/60
trackData = pd.read_csv("track_list.csv")
#trackData = trackData.drop(trackData[(trackData.LEN < lowlim)].index) #remove all cycles below 10 frames length
#trackData = trackData.drop(trackData[(trackData.LEN > lowlim)].index) #remove overlong cycles (for purpose we are interested in)
trackData = trackData.drop(trackData[(trackData.GEN > 2)].index)  #remove generations above 3 as only 1 and 2 are useful in this particular dataset. this will vary between datasets. restrict  to full cell cycle stages (or near full)
#print(trackData)
founder_gen=trackData["GEN"] == 0 #split into generations so we can find family ID's per generation and inspect  branches of family trees sequentially.

print("FOUNDER_GEN!",founder_gen)
gen1=trackData["GEN"] == gen_look_1
gen2=trackData["GEN"] == gen_look_2
#gen3=trackData["GEN"] == gen_look_3

print("founder_gen:",len(trackData[founder_gen]))
founder_gen_len=len(trackData[founder_gen])
gen1_len=len(trackData[gen1])
gen2_len=len(trackData[gen2])
#gen3_len=len(trackData[gen3])

#print(trackData[gen1]["LEN"][1:2].to_string(index=False))
#ref_value=int(trackData[gen1]["LEN"].iloc[2])
#print(int(ref_value)+10)
#a = trackData[:][13]
#print(trackData.iloc[33][1])

fx1list=[] #set up empty lists, these will take tuples containg family ID's,  cell cyclelengths per lineage, generation is encoded in family id by number of  'x's in name. this is the list with the cell cycle length for the founders  generation and the 1st generation
g1x2list=[] #list for 1st generation and 2nd generation
for track in range(0,founder_gen_len):
	founder_family=str(trackData[founder_gen]["FAM_ID"].iloc[track]) #get name of founder family in this loop. 
	cc0_len=int(trackData[founder_gen]["LEN"].iloc[track]) #founder cell cell cycle length
	print("family: {}, cell cycle length: {} ".format(founder_family,cc0_len))
	fam_subdata1 =(trackData[gen1])[(trackData[gen1])['FAM_ID'].str.contains(str(founder_family+"x"))] #look for daughter cells from founder family and subset into dataframe
	print(fam_subdata1)
	for subtrack1 in range(0,len(fam_subdata1)): #iterate through daughter cell  dataframe to extract each daughter cell cell cycle lengths
		gen1_fam=str(fam_subdata1["FAM_ID"].iloc[subtrack1])
		cc1_len=fam_subdata1["LEN"].iloc[subtrack1]
		print("cc1_len",cc1_len)
		genFx1 = (gen1_fam, cc0_len, cc1_len) #ID/comparison tuple: family ID, founder length, gen1 length
		print("genFx1", genFx1)
		fx1list.append(genFx1) #add ID/comparison  tuple to list, to plot as graph
		fam_subdata2 =(trackData[gen2])[(trackData[gen2])['FAM_ID'].str.contains(gen1_fam)]
		print(fam_subdata2)
		for subtrack2 in range(0,len(fam_subdata2)): #repeat steps for founder/gen1 for gen1/gen2
			gen2_fam=str(fam_subdata2["FAM_ID"].iloc[subtrack2])
			cc2_len=fam_subdata2["LEN"].iloc[subtrack2]
			print("cc2_len",cc2_len)
			genFx2 = (gen2_fam, cc1_len, cc2_len)
			print("genFx2", genFx2)
			g1x2list.append(genFx2)

			#gen2_fam=str(fam_subdata2["FAM_ID"].iloc[subtrack2])
			#fam_subdata3 =(trackData[gen3])[(trackData[gen3])['FAM_ID'].str.contains(gen2_fam)]
			#print(fam_subdata3)
	


'''
print(fx1list) #plot founder/gen1 comparison
for a in fx1list:
	if "x3" in a[0]:
		continue
	print(a)
	plt.scatter(a[1],a[2],c="black")
plt.show()
'''	
x_coord_list=[]
y_coord_list=[]
print(g1x2list) #plot gen1/gen2 comparison. :::: need to create output to display both on same multi plot
for a in g1x2list:
	if "x3" in a[0]:
		continue
	print(a)
	if a[1] < lowlimX or a[1] > highlimX or a[2] < lowlimY or a[2] > highlimY:
		continue
	plt.scatter(a[1]*conversionFrameR,a[2]*conversionFrameR,c='black',alpha=0.1)
	x_coord_list.append(a[1]*conversionFrameR)
	y_coord_list.append(a[2]*conversionFrameR)
b,m=polyfit(np.array(x_coord_list),np.array(y_coord_list),1)
	
#plt.plot(np.array(x_coord_list),b+m*np.array(x_coord_list),'-')
plt.title("Correlation between generation {} and generation {} from dataset {}".format(gen_look_1+1,gen_look_2+1,dataset))	
plt.xlabel("generation {} time (mins)".format(gen_look_1+1))
plt.ylabel("generation {} time (mins)".format(gen_look_2+1))
plt.show()


	#gen1Data=trackData[gen1]
	#fam_subdata =gen1Data[gen1Data['FAM_ID'].str.contains(founder_family)]
	#print(fam_subdata)

#below: started setting up drawing a lineage tree. 
'''
image_dimensions=(1000,1000) 
i2 = Image.new("1",image_dimensions,"white") #initialise canvas
mapping = ImageDraw.Draw(i2)
print(image_dimensions[1])
#rectangle_params #describe rectangle parameters. drawing rectangles for 'bars' of cell cycle length. rather than doing bar chart, just designate start and end points of boxes and draw them
start_lineage_x=10
start_lineage_y=((image_dimensions[1])/2)
rect_height=20
gen_gap_1=rect_height*8

mapping.rectangle([(start_lineage_x,start_lineage_y),(start_lineage_x+10*ref_value,start_lineage_y+rect_height)],fill="black")

i2.show()
'''