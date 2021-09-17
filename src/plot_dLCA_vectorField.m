%function to plot LCA displacement field
%OM @ MISL, om82@drexel.edu 20 Nov 2014


function q = plot_dLCA_vectorField(C,d,scl,lSpec)
q = quiver(C(:,1),C(:,2),scl*d(:,1),scl*d(:,2),lSpec);
set(q,'autoScale','off','LineWidth',2)