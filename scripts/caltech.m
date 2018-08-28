clear all;
% ccount1 = 1;
% ccount2 = 1;
quote = '/';
if ispc
    quote = '\';

display('Input: 1. File containing feature vector in each line. 2. File to classify \n Option: append Fisher Vector to given feature vector.');
display('Remember, we need VLFEAT to run the FIsher vector implementation.')

BPTrain = [];
BPTest = [];
trainLabel = [];
testLabel = [];
noline = 0;
filename = uigetfile('*.csv','Enter feature vector file-name: ');
YN = input('Do you want to append Fisher Vector as feature vectors to your input?[y/n]:','s');
if(YN=='Y' || YN=='y')
% 	imgfolder = input('Enter image foldername','s');
	d1 = uigetdir('Folder containing images');
    dirs = dir(d1);
    dirs(1:2) = [];
    dirs = {dirs.name};
    [~,b]= size(dirs);
    for i=1:b
    jpgfiles = dir(fullfile(strcat(d1,strcat(quote,char(dirs(i)))),'*.jpg'));
    [a,~] = size(jpgfiles);
    for j=1:a
        copyfile(strcat(d1,strcat(quote,strcat(char(dirs(i)),strcat(quote,jpgfiles(j).name)))),d1)
    end
    
    end
    jpgfiles = dir(fullfile(d1,'*.jpg'));
end
if(strcmp(filename(end-3:end),'.csv')~=0)
    display('Input needs to be csv');
end
M = csvread(filename);
[noline,~] = size(M);
fisherencode = 1024;
ccount1 = 1;
ccount2 = 1;
trainLabel = [];
numClusters = 16;
p={};
index = 1;
prevlabel = 1;
labelling = 0;
labelTrImg = {};
labelTeImg = {};
encoding = [];
% counteachset = zeros(1,5);
c = 1;
for i=1:noline
%     counteachset(1,c)  = counteachset(1,c) + 1;
if(YN=='Y'|| YN=='y')
    jpgfilenow = jpgfiles(i);
%     xmlname = jpgfilenow.name;
    namm = strcat(d1,strcat(quote,jpgfilenow.name));
    J = imread(namm);
    if(size(J,3)<3)
        I = single((imread(namm))) ;        
    else
        I = single(rgb2gray(imread(namm))) ;
    end
    [f,d] = vl_sift(I) ;
    v = d';
      if(size(v,1)<16)
            vv = vertcat(v,zeros(16-size(v,1),128));
            vv = vv';
      else
          vv =v';
      end
    vv = double(vv);
    [means, covariances, priors] = vl_gmm(vv, numClusters);
    encoding = vl_fisher(vv, means, covariances, priors)';
end  
	display(M(i));
	B = horzcat(encoding,M(i,1:end-1));
	if(YN=='Y'||YN=='y')   
        B = horzcat(encoding,M(i,1:end-1));
    end
	display(B);
    BPTrain = vertcat(BPTrain,B);
	trainLabel = vertcat(trainLabel,M(i,end));

end
 totalTrain = BPTrain;%(randperm(size(BP,1)),:);
 totalTest = M(1,1:end-1);%(randperm(size(BPTest,1)),:); 
 Train = trainLabel;
 Label = M(1,end);
% Mdl = fitctree(BPTrain',trainLabel);
MdlSVM = fitcecoc(BPTrain,trainLabel);
%  YDT = predict(Mdl,totalTest');
 YSVM = predict(MdlSVM,horzcat(encoding,M(1,1:end-1)));
 display(YSVM);
 if(YN=='Y'|| YN=='y')
     jpgfiles = dir(fullfile(d1,'*.jpg'));
     [a,~] = size(jpgfiles);
    for j=1:a
        delete(strcat(d1,(strcat(quote,jpgfiles(j).name))));
    end
 end