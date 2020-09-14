#GDD Clustering
#
#%--------------------------------------------------------------------------
#%Author:  Emre Güngör, Ahmet Özmen
#% Please Cite as (if you use in you study or work)
#%       Güngör, Emre, and Ahmet Özmen. "Distance and density based 
#%       clustering algorithm using Gaussian kernel." Expert Systems with 
#%       Applications 69(2017): 10-20.
#%
#% Please Note:
#%       Standard error of the mean is used in variables DevD,GPD_CPU,
#%       DPD_CPU, which is sqrt(1/N)*standard deviation. It is used for  
#%       estimation of the standard deviation of a sampling distribution.
#%       Reason for that to correctly estimate relation between samples and
#%       region that samples are resides.
#%       
#%       Also SIP deviation calculations are computed using 
#%       std=sqrt((1/N)*(sum(xi^2)-(sum(xi)^2))/N);
#%       for performance and coding ease. Other deviations can be calculated
#%       this way for performance in loops instead of std() function. 
#%   
#%       Code is not optimized - And can contain errors please help us
#%       if you see any errors by sending e-mail to us. Thanks in advance :)
#%       
#% Contact Info:
#%       Emre Güngör: emregungor84@hotmail.com
#%       Ahmet Özmen: ozmen@sakarya.edu.tr
#%--------------------------------------------------------------------------
#
#
#------------------------------------------------------------------------------
#define libraries
import numpy as np

#------------------------------------------------------------------------------
#GDD clustering algorithm convertion code blocks

#function func_gdd (input, output1, output2) #Pass by reference
def func_gdd(dataset):

	#%%PREPROCESSING

	##Clustering Timer Start
	#timeClusteringStart = time.time();

	dataUnfiltered=dataset;
	N=len(dataUnfiltered[:,1]); #Number of Data Samples

	#Remove duplicate values

	##numpy library unique function           
	data, indexUnique = np.unique(dataUnfiltered, axis=0, return_inverse=True)


			
	#Delete Unneded data        
	del dataUnfiltered; #Delete variable names del(obj) delete values 

	#____________________________
	#Return if N=2 as 2 cluster
	if (len(data[:,1])==2 or (len(data[:,1]) <= len(dataset[:,1])/100) ):
		#Also terminate if unique data less than 1 percent of original data
		clusters={};
		seedPoints=[];
		for i in range(len(data[:,1])):
			clusters[i]=[i];
			seedPoints.append=i;
	   
		#Re-order Cluster Indices According to Original DuplicateDataSet if Any
		for i in range(len(clusters[:,1])):
			ClustMemCount               =len(clusters[i,2]);
			orgIndices                  =[];
			SeedPointSelections         =[];
			for j in range(ClustMemCount):
				orgIndices              =[x for x in np.where((indexUnique == clusters[i][j]))[0] ];
				clusters[i][j]=[];
				for k in range(len(orgIndices[:,1])):
					clusters[i].append  =orgIndices[k];
	               
			seedPointSelections         =[x for x in np.where((indexUnique == seedPoints[i]))[0] ];
			seedPointSize               =len(SeedPointSelections[:,1]);
			seedPointIndex              =1;
			if(seedPointSize>1):
				seedPointIndex          =(round(seedPointSize/2));
			seedPoints[i]=seedPointSelections[seedPointIndex];
		#termination of program   
		quit(); #end if data size       = 2 sample (no need clustering)
	#End of Cluster Datasize check
	#____________________________
	 
	#MINIMUM AS GROUND ZERO 
	#mean of min max interval to eliminate shifting artifacts
	N=len(data);
	D=len(data[0]);
	for j in range(D):
		minX=min(data[:,j]);
		for i in range(N):
			data[i,j]=data[i,j]-minX;
		
	N=len(data); #Number of Data Samples
	D=len(data[0]); #Dimension of Data Space    
	#------------------------------------------------------------------------------

	#------------------------------------------------------------------------------

	## MEAN & DEVIATION
	#Mean for all dimensions
	meanD=np.zeros([D,1]);

	#Deviation for all dimensions
	DevD=np.zeros([D,1]);

	#Mean
	for i in range(N):
		for j in range(D):
			meanD[j]=meanD[j]+data[i,j];

	for j in range(D):
		meanD[j]=(meanD[j]/N);


	#Dimensional Deviance Coeff (DevD)= Deviance * sqrt(1/N)
	for i in range(N):
		for j in range(D):
			DevD[j]=DevD[j]+np.power((data[i,j]-meanD[j]),2);

	for j in range(D):
	   DevD[j]=(1/N)*np.sqrt(DevD[j]); 


	#------------------------------------------------------------------------------

	#------------------------------------------------------------------------------
	## REMOVE UNNECESSARY DIMENSIONS
	EraseDimension=[];

	for j in range(D):
		#--------------------
		if( DevD[j]==0 ): #there is no differentiation in this dimension
		   #Erase Dimension Properties 
		   EraseDimension.append(j);
		   


	#------------------------
	#Erase Unnecessary Dimensions for Clustering--------------EXTREME CASES
	for i in range(len(EraseDimension)):
		data=np.delete(data,[j],axis=1); 
		meanD=np.delete(meanD, 1, i); 
		DevD=np.delete(DevD, 1, i); 


	D=len(data[0]); #Dimension of Data Space

	#------------------------------------------------------------------------------

	##Measure Elapsed Time
	#timePreprocessing = time.time()
	#print("Time for Preprocessing"+ str(timePreprocessing-timeFigureWrite))

	#------------------------------------------------------------------------------
	#%% CLUSTERING 
	 
	#%% STEP 1- STATISTICAL PROPERTY COMPUTATIONS
	 
	dm=np.zeros((N,N));
	gm=np.zeros((N,N));
	ag=np.zeros([N]);
	gpm=np.zeros([N]);
	gpd=np.zeros([N]);
	dpm=np.zeros([N]);
	dpd=np.zeros([N]);
	   
	for i in range(N):
		for j in range(N):
			tempDist=0;
			tempExp=0;
			for k in range(D):
				tempDist=tempDist+np.power((data[i,k]-data[j,k]),2);
				tempExp=tempExp+ (np.power((data[i,k]-data[j,k]),2))/(2*np.power(( np.sqrt((meanD[k]*DevD[k])/(2*np.pi)) ),2));

			tempDist=np.sqrt(tempDist);
			dm[i,j]=tempDist;
			gm[i,j]=np.exp(-1*tempExp);
			ag[i]=ag[i]+gm[i,j];

		gpm[i]=(ag[i]-1)/(N-1);
		dpm[i]=np.sum(dm[i,:])/(N-1);

	 
	for i in range(N):
		for j in range(N):
			gpd[i]= gpd[i]+np.power((gm[i,j]-gpm[i]),2); #Variance before taking square root
			dpd[i]= dpd[i]+np.power((dm[i,j]-dpm[i]),2); #Variance before taking square root
			
		gpd[i] = (1/(N))*np.sqrt(gpd[i]);
		dpd[i] = (1/(N))*np.sqrt(dpd[i]);

	# ---------------------------
	spaceGaussDeviance=np.mean(gpd);
	spaceGaussMean=np.mean(gpm);
	spaceDistanceDeviance=np.mean(dpd);
	spaceAGMean=np.mean(ag);
	# ---------------------------    

	##Measure Elapsed Time
	#timeClustering1Step = time.time()
	#print("Time for Clustering Step 1: "+ str(timeClustering1Step - timePreprocessing))

	##PAUSE for Benchmarking and Debugging process   
	#print("Here Line ..."); 
	#input("Press Enter to continue...")


	#%% STEP 2 - CLUSTERING PROCESS  

	#clusters={};
	clusters=[];
	dataSetIndices=list(range(N)); 
	dataSetIndicesMaxGaussIndex=-1; 
	seedPoints=[];

	#♣SEED POINT SELECTION FOR CLUSTERING PROCESS ON REGION GROWTH
	#    take maximum peak energy of gaussian to pinpoint maximum density in space
	maxGaussianMean = np.max(gpm);  
	maxGaussIndex = np.argmax(gpm);
	dataSetIndicesMaxGaussIndex=maxGaussIndex; #At first its same then removed arry

	maxGaussianVal=maxGaussianMean*N; 
	seedIndex=maxGaussIndex;
	seedAccumulativeGaussian=maxGaussianVal;

	##PAUSE for Benchmarking and Debugging process   
	#print("Here Line ____");  
	#input("Press Enter to continue...")
	   
	##LOOP - AFTER FINDING CLUSTER FIND NEXT CLUSTER THAT IS MOST DENSE
	while(len(dataSetIndices)!=0):
			
		#Cluster Centroids
		seedPoints.append(seedIndex);


		####SEED DISTANCE CALCULATIONS
		#--------------------------------------------------------------------------
		#
		#CLUSTER FIXED COMPONENTS (EQ. 8-9)
		fgdt=(gpm[seedIndex]/(np.pi*np.sqrt(gpd[seedIndex])))-(np.pi*gpd[seedIndex]);
		fdt=abs((dpm[seedIndex])/(np.log(gpm[seedIndex]*N+dpd[seedIndex])*(np.log(gpd[seedIndex]))));
		#
		#DEVIANCE FOR EACH NEIGHBOUR RECURSE NEEDED TO MAINTATIN GRADUALITY
		#--------------------------------------------------------------------------    
	 	#NSL creation and Seed (Cluster Centroid) Assignments
		vlusterIndices=[];
		nsl=[]; #Neighbor Search List
		nsl.append(seedIndex);
	   
		
		clusterIndices=[];
		clusterIndices.append(seedIndex);

		#--------------------------------------------------------------------------
		#Remove Cluster Centroid from Unclustered Sample List
		del dataSetIndices[dataSetIndicesMaxGaussIndex];
		#--------------------------------------------------------------------------
		#SPL Mean and Deviation Variables
		#----------------------------------------------------------------------
		#Variance of Neighborhood variances
		meanBase=0; #Ground Zero Mean
		mg=0;
		vg=0; #Ground Zero Variance
		diffcounter=0;
		#Gaussian CounterPart
		gMeanBase=0;gMg=0;gVg=0;
		#----------------------------------------------------------------------

	#    #PAUSE for Benchmarking and Debugging process    
	#    print("Here Line ---"); 
	#    input("Press Enter to continue...")
		
		#Main Clustering Loop (figure 3)
		while True:
	#        del nsl_new;
			nsl_new=[];
	#        nsl_new=np.empty([]);

			
			if (len(nsl)!=0):
				#Inner Loop of Flow-Chart (figure 3)
				for i in reversed(range(len(nsl))):
					nsp=[];  #Neighbor Sample Points
					
					#SPL (Samples In Proximity List) Calculations______________
	#               #Note: SPL members can be found and calculations can be 
	#               #      done afterwards ....
					pDistMean=0;
					pGaussMean=0;
					cntr=0;
					
					for j in range(N): 
						if(dm[nsl[i],j]<=fdt):
							pDistMean=pDistMean  +dm[nsl[i],j];
							pGaussMean=pGaussMean+gm[nsl[i],j];
							cntr=cntr+1;
						#end
					#end
					if(cntr!=0):
						pDistMean =pDistMean /cntr;
						pGaussMean=pGaussMean/cntr;
					#end  
					
					cntr=0;
					pDistDev=0;
					pGaussDev=0;
					for j in range(N):
						if(dm[nsl[i],j]<=fdt):
							pDistDev=pDistDev+np.power((pDistMean-dm[nsl[i],j]),2);
							pGaussDev=pGaussDev+np.power((pGaussMean-gm[nsl[i],j]),2);
							cntr=cntr+1;
						#end     
					#end
					if(cntr!=0):
						pDistDev=np.sqrt(pDistDev/cntr);
						pGaussDev=np.sqrt(pGaussDev/cntr);
						#Variance of Neighborhood Variances
						diffcounter=diffcounter+1;
						meanBase=meanBase+pDistDev;
						vg=vg+np.power(pDistDev,2);
	#                            
						gMeanBase=gMeanBase+pGaussDev;
						gVg=gVg+np.power(pGaussDev,2);
					#end
	#                        
	#               #Std Deviance Using N loop
	#               #-----------------------------
	#%              #std=sqrt((1/N)*(sum(xi^2)-(sum(xi)^2))/N);
	#               #----------------------------
	#               #Calculations below are standard deviations of
	#               #variances
					mg= np.power((meanBase),2)/diffcounter;
					stdVariance=np.sqrt((1/(diffcounter))*np.abs(vg-mg));
					meanVariance=meanBase/diffcounter;
	#                           
					gMg=np.power((gMeanBase),2)/diffcounter;
					gStdVariance=np.sqrt((1/(diffcounter))*np.abs(gVg-gMg));
					
	#                print("End Range N");
					#END for Range N   
				
	#_______    #END for NSL loop
	#______________________________________________________________________________________________________________
	#               #GRADIENT CHANGE EQUATIONS (EQ. 10 -11 )  
					gdt=(meanVariance*(ag[nsl[i]]/ag[seedIndex]))+2*np.pi*stdVariance*(spaceGaussDeviance/gpd[nsl[i]]);
					ggdt=gStdVariance;
	#              
	#               #THRESHOLDS
					thresholdA=fgdt-ggdt; #FGDT-GGDT;
					thresholdB=fdt+gdt;   #FDT+GDT;
					
	#               #___________________________________________________________________________________________
					for k in reversed(range(len(dataSetIndices))):
						if ( gm[nsl[i],dataSetIndices[k]]>=thresholdA and dataSetIndices[k]!=nsl[i] ):
							if( dm[nsl[i],dataSetIndices[k]]<=thresholdB and dataSetIndices[k]!=nsl[i] ):
								nsp.append(dataSetIndices[k]); #point k is in the neighbourhood and same cluster
								del dataSetIndices[k]; #Remove element from list at k
	#
	#               #%ADD Neighbors to Cluster
					if (len(nsp)!=0):
						for x in nsp:
							clusterIndices.append(x);
						for x in nsp:
							nsl_new.append(x);
	#            print("End 1st FOR (NSL) Loop");            
	#           #end LOOP 1st FOR
				del nsl;
				nsl=nsl_new; #New Found Neighbors are now NSL
				
			else:
				#There is no neighbours that fit similarity measures
	#            print("size nsp is 0");
				break;
			
			#print("End 1st IF");       
			#end IF 1st
		   
	#        #PAUSE for Benchmarking and Debugging process    
	#        print("Here Line ---"); 
	#        input("Press Enter to continue...")
		 
	#    print("End While(True)");
		#end-indent While(1)
			 
		#FOR NEXT CLUSTERING - SETTING
		clusters.append(clusterIndices);

		#Set New Maximum Gauss Point - NEXT MOST DENSE POINT  
		maxGaussianMean=0;
		maxGaussIndex=0;
		for i in range(len(dataSetIndices)):
			if(gpm[dataSetIndices[i]]>maxGaussianMean):
				maxGaussianMean=gpm[dataSetIndices[i]];
				maxGaussIndex=dataSetIndices[i];
				dataSetIndicesMaxGaussIndex=i; #Keep track which indices are set
			#end
		#end
		seedIndex=maxGaussIndex;
		
		#Special Case - No Centroid found & 1-Sample Remains
		if(seedIndex==0 and len(dataSetIndices)==1):
			#OPTION: Can add 1 sample to new cluster
			clusters.append([dataSetIndices[0]]);
			seedPoints.append(dataSetIndices[0]);
			del dataSetIndices[0]; #dataSetIndices[1]=[];
			break; 
		#end
		
		#Special Case - No Centroid found & N-Sample Remains
		if(seedIndex==0 and len(dataSetIndices)>1):
			for n in range(len(dataSetIndices)):
				clusters.append(dataSetIndices[0]); #CLUSTERS{end+1}=DataSetIndices(1); 
				seedPoints.append(dataSetIndices[0]);
				del dataSetIndices[0]; 
			#end
		#end
		if(seedIndex==0):
			break; 
	#    #PAUSE for Benchmarking and Debugging process    
	#    print("Here ---"); 
	#    input("Press Enter to continue...")    
	#    print("End While(1)");        
		#end while(1)
	#print("End Whole Clustering");       
	#end -statement #END CLUSTERING

	#    
	##Measure Elapsed Time
	#timeClustering2Step = time.time()
	#print("Time for Clustering Step 2: "+ str(timeClustering2Step-timeClustering1Step))
	
	
	#%% POST-PROCESS STEP
	clustersnew=[[]];
	clustMemCount=[];
	seedPointSize=[];
	seedPointIndex=[];
	#
	for i in range(len(clusters)):
		clustMemCount=len(clusters[i]);
		
		if(len(clustersnew)==i):
	#        clustersnew.append([None]);
			clustersnew.append([]);
		
		
	#%     disp([ 'Cluster Number: ' num2str(i) ' ClusterMember: ' num2str(ClustMemCount)]);
	#%     disp([ 'Total Index Size '  num2str(size(indexC)) ]);
		for j in range(clustMemCount):
			#orgIndices=find(indexUnique(:)== clusters[i][j]); 
			#orgIndices=[x for x in indexUnique if (indexUnique == clusters[i][j])];
			orgIndices=[x for x in np.where((indexUnique ==clusters[i][j]))[0] ];
			
			for k in range(len(orgIndices)):
				#CLUSTERSnew{i}(end+1)=OrgIndices(k);
				clustersnew[i].append(orgIndices[k]);

		seedPointSelections=[x for x in np.where((indexUnique == seedPoints[i]))[0] ];
		
		seedPointSize=len(seedPointSelections);
		seedPointIndex=0;
		if(seedPointSize>1):
			seedPointIndex=(np.round(seedPointSize/2));
		#end
		seedPoints[i]=seedPointSelections[seedPointIndex];
	#end
	clusters=clustersnew;    

	##Measure Elapsed Time
	#timeClusteringEnd = time.time()
	#print("Total Time for Clustering : "+ str(timeClusteringEnd-timeClusteringStart))
	#------------------------------------------------------------------------------

	return clusters,seedPoints;	
#------------------------------------------------------------------------------

