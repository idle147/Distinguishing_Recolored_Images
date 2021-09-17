function [C,d,p] = lca_corners_localDisp_globalParams(I,opts)
nCorners = opts.nCorners;
hW = opts.hW;
delta = opts.delta;
u = opts.u;
cThresh = opts.cThresh;
iRC = opts.iRC;
dGamma = opts.dGamma;
maxiter = opts.maxiter;
convergedThresh = opts.convergedThresh;
%opts.makePlots,opts.displayText

if ~isfield(opts,'NComparison'); %is number of comparison channels
    opts.NComparison = 2;
end

cBlkSize = floor(sqrt((size(I,2)*size(I,1))/nCorners));

mask = imedgemask(I,hW+delta+1); %don't get corners near edges
C = selectiveCorners_withMask(rgb2gray(I),mask,'MinimumEigenValue',cBlkSize,cThresh,nCorners);

if isempty(C)
    switch opts.NComparison
        case 1
            d = [];
            p = [];
        case 2
            d = cell(1,2);
            p = cell(1,2);
    end
    return
end

wbar = waitbar(0,'Estimating local LCA...');
switch opts.NComparison
    case 1
        d = zeros(size(C));
        for ci = 1:size(C,1);
            d(ci,:) = localLCAdisplacement_fixedDS(I,iRC,opts.iCC,C(ci,:),hW,u,delta);
            waitbar(ci/size(C,1),wbar);
        end
        close(wbar);
        p = estimateJFparamsGaussNewton(C,d,size(I),dGamma,maxiter,convergedThresh,0,0);
    case 2
        
        d = cell(1,2);
        for ci = 1:size(C,1);
            [d{1}(ci,:),d{2}(ci,:)] = localLCAdisplacement_fixedDS_2Comp(I,iRC,C(ci,:),hW,u,delta);
            waitbar(ci/size(C,1),wbar);
        end
        close(wbar);
        p{1} = estimateJFparamsGaussNewton(C,d{1},size(I),dGamma,maxiter,convergedThresh,0,0);
        p{2} = estimateJFparamsGaussNewton(C,d{2},size(I),dGamma,maxiter,convergedThresh,0,0);
end