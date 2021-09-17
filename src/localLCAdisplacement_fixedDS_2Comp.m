%Script to estimate local LCA displacement at an image location
%Finds (fractional) pixelwise shift of reference channel that results in
%maximum similarity to comparsison color channel.
%Uses a diamond search method ala source from below

%OM @ MISL, om82@drexel.edu 26 Feb 2015
%Revision 1, 19 May 2015, OM - measure two comparison channels
%simultaneously to save computation required during upsampling process

%Adapted from "Efficient estimation and large-scale evaluation of lateral
% chromatic aberration for digital image forensics" By Gloe et al. 2010
%AND "A New Diamond Search Algorithm for Fast Block-Matching Motion
%Estimation" by Zhu & Ma 2000

%INPUTS:
%I is input image
%iRC, reference channel index (1, 2 or 3 for R, G or B);
%C, corner pt [i,j]
%hW, window half width
%u, upsample factor, search resolution (searches in 1/u steps)
%delta, maximum allowed displacement

%OUTPUTS
%Dhat   - local LCA displacement estimate
%niter  - number of similarity calcs used


function [Dhat1,Dhat2, iCC, niter] = localLCAdisplacement_fixedDS_2Comp(I,iRC,C,hW,u,delta)

%% SETUP
%refernce channel window indices;
%hW*2 square window about a corner.
%+delta on either side to allow for sliding
inds_RefWin = C(1) + ((-1*hW-delta):(hW+delta-1)) ;
jnds_RefWin = C(2) + ((-1*hW-delta):(hW+delta-1));
%Get windowed data from image
Iref_win = I(jnds_RefWin,inds_RefWin,iRC);

iCC = setdiff(1:3,iRC);

%comparison channel window indices
%hW*2 square window about a corner.
inds_CompWin = C(1) + ((-1*hW):(hW-1));
jnds_CompWin = C(2) + ((-1*hW):(hW-1));


%Upsample windowed images
%imresize defaults to bicubic interpolation
Iref_up = imresize(Iref_win,u);

%define shifts in pixel space
delta2 = -1*delta:(1/u):delta; %set of posible displacements
Ndelta2 = length(delta2); %number of candidate shifts in each dim

% DEFINE SEARCH PATTERNS
LDSP = [...
    0, -2;
    1, -1;
    2,  0;
    1,  1;
    0,  2;
    -1,  1;
    -2,  0;
    -1, -1];

SDSP = [...
    0, -1;
    1,  0;
    0,  1;
    -1,  0];

niter = [0,0];
for iComp = 1:2;
    %Get windowed data from image
    Icomp_win = I(jnds_CompWin,inds_CompWin,iCC(iComp));
    %upsample
    Icomp_up = imresize(Icomp_win,u);
    
    
    
    %% INITIALIZE SEARCH COMP LAYER 1
    di = delta*u+1; %start at middle
    dj = delta*u+1;
    
    %get indicces of Ref block
    [indsRef, jndsRef] = getIndsRef(di,dj,hW,u);
    
    dmax = corr2(Iref_up(jndsRef,indsRef),Icomp_up); %measure correlation
    niter(iComp) = 1;
    
    % d = [di, dj];
    ds = [];
    
    %% Large Diamond Search pattern
    
    isLDSP = 1; %keep doing this until we converge upon a particular area
    while isLDSP;
        
        %make sure no overlap between previous LDSP
        if niter(iComp) == 1; %first time? no overlap
            d = ones(length(LDSP),1)*[di,dj]+LDSP;
        else %otherwise, probably some overlap
            d = setdiff(ones(length(LDSP),1)*[di,dj]+LDSP,ds,'rows');
        end
        
        %make sure no points exceed max/min
        %remember these are indices
        d(d > (2*u*delta+1)) = (2*u*delta+1);
        d(d < 1) = 1;
        
        %initilaze similarity vec we are about to calculate
        C1 = zeros(1,size(d,1));
        
        %iterate over all new displacements
        for i1 = 1:size(d,1)
            %measure correlation
            C1(i1) = corr2(Iref_up(jndsRef+d(i1,2)-dj,indsRef+d(i1,1)-di),Icomp_up);
            niter(iComp) = niter(iComp)+1; %keep track of number of correlations
        end
        
        if max(C1) > dmax %does this LDSP yield a new max?
            [dmax, idmax] = max(C1); %save the max
            ds = [d;ds]; %keep track of used displacemements
            di = d(idmax,1); %new center
            dj = d(idmax,2); %new center
            %get indices of Ref block at new center
            [indsRef, jndsRef] = getIndsRef(di,dj,hW,u);
        else %if LDSP does not yield a new max
            isLDSP = 0; %switch to SDSP
        end
        
    end
    
    %%
    C2 = zeros(1,length(SDSP));
    
    ds2 = [SDSP(:,1)+di SDSP(:,2)+dj]; %search pts
    
    d = ones(length(SDSP),1)*[di,dj]+SDSP;
    
    %prevent from going over max/min search
    d(ds2 > (2*u*delta+1)) = (2*u*delta+1);
    d(ds2 < 1) = 1;
    
    for i2 = 1:length(SDSP)
        C2(i2) = corr2(Iref_up(jndsRef+d(i2,2)-dj,indsRef+d(i2,1)-di),Icomp_up); %measure correlation
        niter(iComp) = niter(iComp)+1;
    end
    if max(C2) > dmax
        [~, idmax] = max(C2);
        di = d(idmax,1);
        dj = d(idmax,2);
    end
    
    switch iComp
        case 1
            Dhat1 = delta2([di, dj]);
        case 2
            Dhat2 = delta2([di, dj]);
    end
end
end

function [indsRef, jndsRef] = getIndsRef(di,dj,hW,u)
indsRef = (di-1) + (1:2*hW*u); %x indices of Upsampled Ref image centered at di
jndsRef = (dj-1) + (1:2*hW*u); %y indices of Upsampled Ref image centered at dj
end