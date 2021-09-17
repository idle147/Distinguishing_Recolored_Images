function p = lca_paramJFfrompair(C1,C2,dhat1,dhat2)

%build Matrix
A = [1 0 dhat1(1);
    0 1 dhat1(2);
    1 0 dhat2(1)];
% keyboard;
%build vector son
b = [C1(1); C1(2); C2(1)];

x = A\b; %solve for LCA parameters

p = [x(1:2); 1-(1/x(3))]; %put into parameter form