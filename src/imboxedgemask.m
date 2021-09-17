%BOX MASK
function boxMask = imboxedgemask(I,x,y,w,h,hT) 

% hT = 5; %pixel boarder to not look at corners (half thickness)
boxMask = ones(size(I,1),size(I,2));

%top edge
boxMask(y+(-hT:hT),x+(-hT:w+hT)) = 0;
%bottom edge
boxMask(y+h+(-hT:hT),x+(-hT:w+hT)) = 0;
%left edge
boxMask(y+(0:h),x+(-hT:hT)) = 0;
%bottom edge
boxMask(y+(0:h),x+w+(-hT:hT)) = 0;

boxMask = logical(boxMask);