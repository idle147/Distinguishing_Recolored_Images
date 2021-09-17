%PTSINBOX function that returns indices of points that are in a specified box

%OM @ MISL, om82@drexel.edu, 30Jan2015

%INPUTS:
%   C   -   Coordinates of points to test, first col is x second col is y
%       (size N x 2)
%   Box -   Dimensions of the box [x y w h], x y specify bottom left corner
%   IncludeEdges -  logical that indicates whether to include the edges as
%       in the box (1) or disclude the edges as part of the box (0)

%OUTPUTS:
%   iCInBox     - indices of points that are in the box  
%   iCNotInBox  - indices of points that are not in the box
%   isBox       - N x 1 vector of logicals. 1 indicates in box, 0 not.

function [iCInBox,iCNotInBox,isBox] = PtsInBox(C,Box,IncludeEdges)

if nargin < 3
    IncludeEdges = 0; %don't include edges on default
end

if Box(3) < 0 || Box(4) < 0; %make sure dimensions are good
    error('PTSINBOX: Box width/height must be positive')
end


%% FILTER PTS
if ~IncludeEdges %Don't Include Edges
    i1 = find(C(:,1) > Box(1)); %pts greater than left edge
    i2 = find(C(:,1) < (Box(1)+Box(3))); %pts smaller than right edge
    i3 = find(C(:,2) > Box(2)); %pts greater than bottom edge
    i4 = find(C(:,2) < (Box(2)+Box(4))); %pts less than top edge
else %Include Edges
    i1 = find(C(:,1) >= Box(1));
    i2 = find(C(:,1) <= (Box(1)+Box(3)));
    i3 = find(C(:,2) >= Box(2));
    i4 = find(C(:,2) <= (Box(2)+Box(4)));
end

%% SORT BY IN BOX OR NOT
iCInBox = intersect(i1,intersect(i2,intersect(i3,i4))); %find intersection
iCNotInBox = setdiff(1:size(C,1),iCInBox)'; %find compliment

%% MAKE LOGICAL VECTOR
isBox = zeros(length(C),1);
isBox(iCInBox) = 1;
isBox = logical(isBox); %convert to logical
