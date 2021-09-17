%function to make boxes
% bw = 320;
% bh = 320;
% xOvrlp = 160;
% yOvrlp = 160;

function [INDS,JNDS] = imoverlappingboxes(I,bw,bh,xOvrlp,yOvrlp)

%we can make this prettier/more elegant
ibox = 1;
for ii = 1:(bw-xOvrlp):size(I,2)-bw+1; %starting index of box
    for jj = 1:(bh-yOvrlp):size(I,1)-bh+1; %starting jndex of box
        INDS(ibox,:) = ii + (0:bw-1); %indices of this box
        JNDS(ibox,:) = jj + (0:bh-1); %jndices of this box
        ibox = ibox+1; %increment the box
    end
end