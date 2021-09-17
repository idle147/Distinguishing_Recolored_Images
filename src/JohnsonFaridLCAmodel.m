%function to generate point locations in a comparison channel (and  
%associated displacements) using Johnson and Farid's LCA model

%OM @ MISL, om82@drexel.edu, 18 Nov 2014

%Adapted from:
%Johnson & Farid "Exposing Dig. Forgeries Through Chrom. Aberration" (2006) 
%Gloe et al. "Efficient Est. ... of LCA for Digital Image Forensics" (2010)

%INPUTS:
%R is image coordinate in reference channel @ which to model LCA
%   size is Nx2, where (n,1) is nth coordinate in x (or i), and 
%   (n,2) is nth coordinate in y (or j)
%p is Johnson and Farid LCA parameters, size is 1x3
%   p(1) is optical center in x (or i)
%   p(2) is optical center in y (or j)
%   p(3) is expansion/compression coef, alpha

%OUTPUTS:
%C is image coordinate of comparison channel, corresponding to points
%in input R. size is Nx2
%D is displacement vector. size Nx2

function [C, D] = JohnsonFaridLCAmodel(R,p)

%Calculate position in comparison channel
%EQ (1) in Gloe et al. EQ (4,5) in J&F
C(:,1) = (p(3) .* (R(:,1)-p(1))) + p(1); %x component
C(:,2) = (p(3) .* (R(:,2)-p(2))) + p(2); %y component

%Calculate modeled displacement
%EQ (3,4) in Gloe et. al.
D = R - C;


