# -*- coding: utf-8 -*-
"""
Created on Thu May  3 12:12:32 2018

@author: Emre Güngör
"""



#------------------------------------------------------------------------------
#Library & Dependencies
import scipy.io as sio
import numpy as np

#matlab plot library
import matplotlib.pyplot as plt
#import matplotlib as mpl
#from matplotlib import pyplot as plt
#interactive plotting off
#plt.ioff()

#For color
import matplotlib

#Timer
import time
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#DATA INPUT

#Control Timer
timeStart = time.time()
print("Timer Started.")

#Read mat input file
matFilename='SyntheticData.mat';
matFile = sio.loadmat(matFilename);
#Information on mat file
#sio.whosmat('SyntheticData.mat');

#Take dataset lists from mat construct
datasetList=matFile['XX'];

#Data series taken from dataset
datasetNo=4;
dataset=datasetList[0][datasetNo];
datasetSize=len(dataset);

#Delete Unneded data        
del datasetList; #Delete variable names del(obj) delete values 
del matFile;

#Measure Elapsed Time
timeFileLoad = time.time()
print("Time for Loading File: "+ str(timeFileLoad-timeStart))

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Plot DATA
fig1 =plt.figure(figsize=(12, 12));  # figure & size declaration
plt.plot(dataset[:,0],dataset[:,1],'k.'); # draw figure

##Figure texts
plt.title('2D Spatial Dataset Values');
plt.ylabel('Y');
plt.xlabel('X');
plt.grid(True);
plt.axis('equal'); #Equal Axis values
#plt.axes().set_aspect('equal', 'datalim'); #Equal Aspect Ratio
#plt.suptitle('This is a somewhat long figure title', fontsize=16)

##Save Figure
fig1.savefig(matFilename[0:len(matFilename)-4] + '_' +
                         'dataset'+ '_' +str(datasetNo) +'.png',
                         format='png', dpi=fig1.dpi); # This does, too
#plt.show();
plt.close(); #Close figure #Supresses figure shown on console

#Measure Elapsed Time
timeFigureWrite = time.time()
print("Time for Writing Figure File: "+ str(timeFigureWrite-timeFileLoad))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Call GDD function

#output definitions
clusters=[];
seedPoints=[];
from GDD import func_gdd; #Add file and its function inside to script

timeBeforeGDDcluster = time.time()

#Function Call
clusters,seedPoints=func_gdd(dataset);

#Measure Elapsed Time
timeAfterGDDcluster = time.time()
print("Time for GDD Clustering"+ str(timeAfterGDDcluster-timeBeforeGDDcluster))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Display Cluster Results
plotx=dataset[:,0];
ploty=dataset[:,1];
#color=[0]*len(dataset);
colorVals= [0 for i in range(len(dataset))];
plotColor= [0 for i in range(len(dataset))];

#colorIndex=(cm.rainbow(np.linspace(0, 1, len(clusters))));
#colorIndex=['red','green','blue','cyan','purple'];
colorIndex = matplotlib.cm.rainbow(np.linspace(0, 1, len(clusters)));


for i in range(len(clusters)):
    for j in range(len(clusters[i])):
        colorVals[clusters[i][j]]=i;
        plotColor[clusters[i][j]]=colorIndex[i];



fig2 =plt.figure(figsize=(12, 12));  # figure & size declaration
#plt.scatter(plotx, ploty,color);
#plt.plot( plotx,ploty,
#            color=plotColor,
#            x_label = 'property x',
#            y_label = 'property y',
#            title = 'GDD Cluster Result');
##Figure texts
plt.title('GDD Clustering Result');
plt.ylabel('Y');
plt.xlabel('X');
plt.grid(True);
plt.axis('equal'); #Equal Axis values

plt.scatter(plotx, ploty,color=plotColor,);
#plt.plot( plotx,ploty);
#Save Figure
fig2.savefig(matFilename[0:len(matFilename)-4] + '_' +
                         'GDDresult_' +str(datasetNo) +'.png',
                         format='png', dpi=fig1.dpi); # This does, too

#plt.show();
plt.close();

#------------------------------------------------------------------------------





