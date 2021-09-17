%function to find corners
%!TODO: add method to set min distance btwn corners.

%Change Log:
%+10 Nov 2014 Added functionality to limit number of corners
%+5 Dec 2014 added functionality to only select corners in a masked region
%returned
%+Fixed case where no corners exceed specified metric - OM 27feb2015

function [iC, vC] = selectiveCorners_withMask(I,mask,type,maxD,qThresh,nCorners)
CM = cornermetric(I,type); %Array (size(I)) of corner metrics
CM = CM.*mask;
A = zeros(size(CM)); %initialize "best" corner matrix
for jj = 1:maxD:size(A,1)-maxD;
    for ii = 1:maxD:size(A,2)-maxD;
        [vMax, iMax] = max2d(CM(jj:jj+maxD-1,ii:ii+maxD-1)); %find max corner in window
        A(jj+iMax(1)-1,ii+iMax(2)-1) = vMax; %Save the metric of max corner to matrix
    end
end
[B, iB] = sort(reshape(A,1,[]),'descend'); %sort best corner pts by quality
%get location of sorted best corners. Get first N = nCorners that exceed
%metric value Q > qThresh
if ~isempty(iB(find(B > qThresh, nCorners)))
    [iC(:,2), iC(:,1)] = ind2sub(size(A),iB(find(B > qThresh, nCorners)));
    vC = CM(sub2ind(size(A),iC(:,2),iC(:,1))); %corner quality vals at each pt
else
    iC = [];
    vC = [];
end
end

function [maxval, imax] = max2d(A)
%function to find max value and indices of a 2D matrix
%Owen Mayer @ MISL, om82@drexel.edu 4 Nov 2014

%INPUTS:
%A is a 2 dimensional matrix

%OUTPUTS
%maxval is the maximum value of the matrix
%imax is a vector (L = 2) of the indices of the max value

%if there are repeated maximum values, only (the first, columnwise then
%rowise) one is reported.

[maxval, idx] = max(reshape(A,1,[])); %convert to 1d and max
[imax(1),imax(2)] = ind2sub([size(A,1) size(A,2)],idx); %convert indices back to 2d
end


