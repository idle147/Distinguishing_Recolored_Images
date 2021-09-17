%function to measure error btwn local estimate LCA and Johnson Farid Model
%OM @ MISL, om82@drexel.edu, 18 Nov 2014

%Adapted from:
%Gloe et al. "Efficient Est. (...) of LCA for Dig. Image Forensics" (2010)

%INPUTS:
%R is image coordinate in reference channel @ which to calculate LCA error
%   size is Nx2, where (n,1) is nth coordinate in x (or i), and 
%   (n,2) is nth coordinate in y (or j). Where N is number of coordinates

%Dmeasured is measured LCA displacement at each R(n,:)
%   use localLCAdisplacement.m to measure (adapted from Gloe)

%p is Johnson and Farid LCA model parameters, size is 1x3
%   p(1) is optical center in x (or i)
%   p(2) is optical center in y (or j)
%   p(3) is expansion/compression coef, alpha


%OUTPUTS:
%e is error btwn J&F model and measured LCA, size is Nx2
%e(n,1) is pixel displacement difference between model and measured
%comparison channel location in x dimension (Model - measure) at nth
%location
%e(n,2) is similar but in y direction

function e = eLCA_JF(R,Dmeasured,p)

%EQ (4) from Gloe
[~,Dmodel] = JohnsonFaridLCAmodel(R,p);

%EQ (11) from Gloe
e = Dmodel - Dmeasured;
