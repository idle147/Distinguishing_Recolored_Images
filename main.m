%%  MAIN.M
% Feature Extraction
% 1. Modify the path of the dataset
% 2. The result is stored in 'result'

%%  References
% Mayer, Owen, and Matthew C. Stamm. 
% Accurate and Efficient Image Forgery Detection Using Lateral Chromatic Aberration.
% IEEE Transactions on Information Forensics and Security (2018).


addpath('/src/') %Add source code to path
data_path = dir("The file path of the dataset");
[row, col] = size(data_path)
result = zeros(col, 6); 

for name = 1:row
    %% 0 Parameters 
    HW = 16; %Local LCA measurement window half-width (HW)
    DELTA = 5; %Local LCA measurement maximum displacement
    U = 10; %Local LCA measurement Upsample factor (U=10 if often a better choice but U=5 is used here to)
    thres = 0.0002; 
    corner_num = 5000; 
    
    %% 1 Load Image 
    str = data_path(name, 2) + data_path(name, 1);
    I = imread(str);

    %% 2 Select Key Points (Corners)  
    mask = imedgemask(I,HW+DELTA+1); 
    K = selectiveCorners_withMask(rgb2gray(I),mask,'MinimumEigenvalue',2*HW,thres,corner_num);

    %% 3 Measure Local LCA
    %Initialize 2D vectors to hold local LCA displacement measurements
    dhat = zeros(size(K,1),2); %Green->Red LCA
    %wbar = waitbar(0,'Measuring Local LCA'); %progress bar
    
    for ii = 1:size(K,1) %iterate through each keypoint
        dhat(ii,:) = localLCAdisplacement_fixedDS(I,2,1,K(ii,:),HW,U,DELTA); %use diamond search (Sec. III)
        %waitbar(ii/size(K,1),wbar); %update progress bar
    end
    %close(wbar)

    num = 0; 
    locallength = 0; 
    l = 0;
    con = zeros(size(K,1),1);
    for h = 1:size(K,1)
        con(h,1) = abs(dhat(h ,1)-dhat(h ,2));
        if con(h,1) > thres
            locallength = locallength + con(h,1);
            num = num+1;
        end
    end
      
    %% 4 Estimate Global Fit

%     %Find global model parameters (theta) that best fit the local model
      theta = estimateJFparamsGaussNewton(K,dhat,size(I),0.01,1000,10^-5,0,0); 
%     %Use Least Squares Gauss Newton Method
%     %Theta is 3D with elements: x position of optical center, y posistion of
%     %optical center, and expansion coefficietion. Eq (5)
      [~, d] = JohnsonFaridLCAmodel(K,theta); %Global displacement model (Eq (2))
      globallength  = sum(abs(d(:,1)-d(:,2)));
      htemp = histogram(con,0:0.05:1);
      htemp.Normalization = 'probability';
      htemp.BinWidth = 0.05;
      feature = sort(htemp.Values, 'descend');
      hmax = feature(1,1);
      hdiff = hmax - (feature(1,2) + feature(1,3))/2;
      hvar = var(feature); 
     
     %%
     result(name,:) = [num, locallength , globallength, hmax, hdiff, hvar];
     disp(name);
     save result.mat result;
     close all;
     clear
     load result.mat;
end