clear all
clc

ClustNo=[];   %Unknown Cluster from Database - If clusters known overwrite

%% 2D Noise
mu = [2 3];
SIGMA = [1 1.5; 1.5 3];
rng('default');  % For reproducibility
r = mvnrnd(mu,SIGMA,100);
% figure;
% plot(r(:,1),r(:,2),'+')
X=r;

% n = 100;
% X = randn(n,2);
% % Plot the points.
% plot(X(:,1),X(:,2),'.','markerfacecolor','k');

%% PREPROCESSING

%Remove Duplicates 
% [SortedX T  X]=unique(X,'rows');
% % http://www.mathworks.com/matlabcentral/answers/63229-how-to-find-unique-pages-in-a-3d-matrix
[n,m,p]=size(X);
a=reshape(X,n,[],1);
b=reshape(a(:),n*m,[])';
% c=unique(b,'rows','stable')';
c=unique(b,'rows')';
X=reshape(c,n,m,[]);

%% CLUSTERING
[ CLUSTERS, CentroidPoints] = func_GDD( X );

maxClustNo= size(CLUSTERS,2);
N=size(X,1);
D=size(X,2);

if(size(CLUSTERS,2)==0 && size(CentroidPoints,2)==0)
    disp('No cluster found exiting');
    return;
end

%% FIGURES
% _________________________________________________________________________
% D=2;

    %MARKER INDEX
    markerIndex{1}='o';
    markerIndex{2}='+';
    markerIndex{3}='x';
    markerIndex{4}='s';
    markerIndex{5}='d';
    markerIndex{6}='^';
    markerIndex{7}='v';
    markerIndex{8}='>';
    markerIndex{9}='<';
    markerIndex{10}='p';
    markerIndex{11}='h';
    markerIndex{12}='*';
    for i=13:50
        markerIndex{i}='.';
    end

% %Figure Creation of Clusters for 4D
if(D==4)
    ColorSpecArray=zeros(size(CLUSTERS,2),3);
    
    %X-Y-Z
    figure,
    for i=1:size(CLUSTERS,2)
        ColorSpec = [ rand(), rand(), rand() ];
        ColorSpecArray(i,:)=ColorSpec;
        for j=1:size(CLUSTERS{i},2)
            plot3(X(CLUSTERS{i}(j),1),X(CLUSTERS{i}(j),2),X(CLUSTERS{i}(j),3),'Color',ColorSpecArray(i,:) , 'Marker', markerIndex{i});
            hold on
        end
    end

    grid on;

   
end
% _________________________________________________________________________
% %Figure Creation of Clusters for 3D

if(D==3)


figure,

for i=1:size(CLUSTERS,2)
    ColorSpec = [ rand(), rand(), rand() ];
    for j=1:size(CLUSTERS{i},2)
%         plot3(X(CLUSTERS{i}(j),1),X(CLUSTERS{i}(j),2),X(CLUSTERS{i}(j),3),'Color',ColorSpec , 'LineStyle', 'o');
        plot3(X(CLUSTERS{i}(j),1),X(CLUSTERS{i}(j),2),X(CLUSTERS{i}(j),3),'Color',ColorSpec , 'Marker', 'o');
        hold on
    end
end
   

grid on;
hold off
end
% _________________________________________________________________________
% %Figure Creation of Clusters for 2D
if(D==2)
figure,
axes('position',[.05 .05 .92 .92]);
for i=1:size(CLUSTERS,2)
    ColorSpec = [ rand(), rand(), rand() ];
    for j=1:size(CLUSTERS{i},2)
%         plot(X(CLUSTERS{i}(j),1),X(CLUSTERS{i}(j),2),'Color',ColorSpec , 'LineStyle', 'o');
%         plot(X(CLUSTERS{i}(j),1),X(CLUSTERS{i}(j),2),'Color',ColorSpec , 'Marker', 'o');
        plot(X(CLUSTERS{i}(j),1),X(CLUSTERS{i}(j),2),'Color',ColorSpec , 'Marker', markerIndex{mod(i,12)+1});
        hold on
    end
end
    
for i=1:size(CentroidPoints,2)
    plot(X(CentroidPoints(i),1),X(CentroidPoints(i),2),'*r');
    hold on
end

grid on;
hold off
end
% _________________________________________________________________________
% 1D Figures
if(D==1)
figure,

for i=1:size(CLUSTERS,2)
    ColorSpec = [ rand(), rand(), rand() ];
    for j=1:size(CLUSTERS{i},2)
%         plot(X(CLUSTERS{i}(j),1),0,'Color',ColorSpec , 'LineStyle', 'o');
        plot(X(CLUSTERS{i}(j),1),0,'Color',ColorSpec , 'Marker', 'o');
        hold on
    end
end
    
for i=1:size(CentroidPoints,2)
    plot(X(CentroidPoints(i),1),0,'*r');
    hold on
end

grid on;
hold off
end

% _________________________________________________________________________
% Figure Creation of Clusters for MORE THAN 5D using 3D figures
if(D>5)

figure,

for i=1:size(CLUSTERS,2)
    ColorSpec = [ rand(), rand(), rand() ];
    for j=1:size(CLUSTERS{i},2)
        plot3(X(CLUSTERS{i}(j),1),X(CLUSTERS{i}(j),2),X(CLUSTERS{i}(j),3) ,'Color',ColorSpec , 'Marker', 'o');
%         plot3(X(CLUSTERS{i}(j),1),X(CLUSTERS{i}(j),2),X(CLUSTERS{i}(j),3) ,'Color',ColorSpec , 'LineStyle', 'o');
%         plot3(X(CLUSTERS{i}(j),2),X(CLUSTERS{i}(j),3),X(CLUSTERS{i}(j),9),'Color',ColorSpec , 'LineStyle', 'o');
        hold on
    end
end
   

grid on;
hold off
end

%Correctly see distribution of data points
axis equal
%
% axis equal
% _________________________________________________


