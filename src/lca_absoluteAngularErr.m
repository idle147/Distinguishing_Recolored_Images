%LCA_ABSOLUTEANGULARERR - function to calculate the absoluate angular error between 
%local and global lca displacement vectors.
%if local lca displacement has zero magntidue, a NaN is returned

%INPUTS:
%dglob - Nx2 vector of global lca displacement estimates. (N observations)
%        columns are and x and y components, respectively
%dloc - Nx2 vector of local lca displacement estimates

%OUTPUTS:
%absAng - Nx1 vector of absolute angular errors. Any NaNs occur when dloc
%has zero magnitude - and thus an angular difference can't be measured.

%requires function diffAng - in MISL LCA toolset

%OM@MISL - om82@drexel.edu - 13Apr2015

function absAng = lca_absoluteAngularErr(dglob,dloc)

%absolute angular error between global and local estimates
absAng = abs(diffAng(dglob,dloc)); 

%if local (or global, though unlikely) estimates are zero magnitude, set to
%NaN
absAng(dloc(:,1) == 0 & dloc(:,2) == 0) = NaN;
absAng(dglob(:,1) == 0 & dglob(:,2) == 0) = NaN;