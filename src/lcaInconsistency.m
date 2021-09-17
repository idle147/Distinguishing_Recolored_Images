function E = lcaInconsistency(C,dhat_r,p_r,dhat_b,p_b)

%GLOBAL MODEL DISPLACEMENTS
[~, d_r] = JohnsonFaridLCAmodel(C,p_r);
[~, d_b] = JohnsonFaridLCAmodel(C,p_b);

%ERROR TERM
e1 = (d_r-dhat_r)/p_r(3);
e2 = (d_b-dhat_b)/p_b(3);

E = horzcat(e1,e2); %matricize
