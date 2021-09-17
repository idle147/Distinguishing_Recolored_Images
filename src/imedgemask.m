%create a mask around edge of image, with thickness T
%returns mask only 

function edgeMask = imedgemask(I,T)

edgeMask = ones(size(I,1),size(I,2)); %set mask to one

%set edges to zero
edgeMask(1:T,:) = 0; %top
edgeMask(end-T:end,:) = 0; %bottom
edgeMask(:,1:T) = 0; %left
edgeMask(:,end-T:end) = 0; %right

edgeMask = logical(edgeMask); %make logical