from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd

intervalTime=350
framesNo=65


dataframe=pd.read_csv("appendcsv.csv")
XgroupbyTP=np.array(dataframe.groupby("T")["X"].apply(list))
YgroupbyTP=np.array(dataframe.groupby("T")["Y"].apply(list))
ZgroupbyTP=np.array(dataframe.groupby("T")["Z"].apply(list))

def animate(frame,X,Y,Z):
	ax.clear()
	ax.set_xlim3d([0,750])
	ax.set_xlabel("X (microns)")
	ax.set_ylim3d([0,750])
	ax.set_ylabel("Y (microns)")
	ax.set_zlim3d([0,300])
	ax.set_zlabel("Z (microns)")
	title=ax.set_title("Position of feature detection hits over time, in zebrafish embryo")
	ax.plot(X[frame],Y[frame],Z[frame],linestyle='', marker="o")

fig=plt.figure()
ax=fig.add_subplot(projection="3d")

ani = FuncAnimation(fig,animate,frames=int(framesNo),fargs=(XgroupbyTP,YgroupbyTP,ZgroupbyTP),interval=int(intervalTime))
plt.show()
plt.clf()
