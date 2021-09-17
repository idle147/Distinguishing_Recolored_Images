%Script to measure local LCA displacement at a corner pt
%Finds (fractional) pixelwise shift of reference channel that results in 
%maximum similarity to comparsison color channel.

%OM @ MISL, om82@drexel.edu 07 Nov 2014

%Adapted from "Efficient estimation and large-scale evaluation of lateral
% chromatic aberration for digital image forensics" By Gloe et al. 2010

%INPUTS:
%I is input image
%iRC, reference channel index (1, 2 or 3 for R, G or B);
%iCC, comparison channel index
%C, corner pt [i,j]
%hW, window half width
%u, upsample factor, search resolution (searches in 1/u steps)
%delta, additional search space in pixels
%diagnositics, 0 or 1 for diagnostic plots (does nothing atm - OM 7Nov14)

%OUTPUTS
%CC, correlation matrix at each candidate shift, size
%length(delta2) x length(delta2)
%delta2, candidate shift points

%!TODO: Maybe implement other similarity metrics other than correlation.

function [d, CC, delta2] = localLCAdisplacement_fullsearch(I,iRC,iCC,C,hW,u,delta,diagnostics)

%%%DUMMY TEST parameters. Ignore 
% % close all; clear all;
% % load('subPxLCAtestImage_06-Nov-2014.mat'); %load image
% % I = uint8(I);
% % iRC = 2; iCC = 1; hW = 32; u = 5; delta = 3;
% % %find corners
% % C = corner(I(:,:,2),'MinimumEigenValue'); %corner method on green channel (unshifted)
% % figure; imshow(I);hold on; %show image
% % plot(C(:,1),C(:,2),'co','LineWidth',2) %show corners in cyan
% % 
% % A = C(1,:); clear C; C = A; clear A;
%% Get windowed images

%refernce channel window indices;
%hW*2 square window about a corner. 
%+delta on either side to allow for sliding
inds_RefWin = C(1) + ((-1*hW-delta):(hW+delta-1)) ;
jnds_RefWin = C(2) + ((-1*hW-delta):(hW+delta-1));

%comparison channel window indices
%hW*2 square window about a corner. 
inds_CompWin = C(1) + ((-1*hW):(hW-1));
jnds_CompWin = C(2) + ((-1*hW):(hW-1));

%Get windowed data from image
Iref_win = I(jnds_RefWin,inds_RefWin,iRC);
Icomp_win = I(jnds_CompWin,inds_CompWin,iCC);

%% Upsample windowed images
%imresize defaults to bicubic interpolation
Iref_up = imresize(Iref_win,u);
Icomp_up = imresize(Icomp_win,u);

%% define shifts in pixel space
delta2 = -1*delta:(1/u):delta;
Ndelta2 = length(delta2); %number of candidate shifts in each dim 

%% measure correlation
CC = zeros(Ndelta2,Ndelta2); %preallocate corr (similarity) matrix
for dj = 1:Ndelta2; %iterate through all shift candidates in y
    for di = 1:Ndelta2; %iterate through all shift candidates in x
        indsRef = (di-1) + (1:2*hW*u); %x indices of Upsampled Ref image for candidate shift
        jndsRef = (dj-1) + (1:2*hW*u); %y indices of Upsampled Ref image for candidate shift
        CC(di,dj) = corr2(Iref_up(jndsRef,indsRef),Icomp_up); %measure correlation
    end
end

switch diagnostics
    case true
        %space for some diagnistic figure generation etc
        h = figure; close(h); %dummy figure for now
    case false
end


%% Calculate offset

%find max similarity shift
[~,imax] = max2d(CC);

%convert to pixel dimensions
d = delta2(imax);
