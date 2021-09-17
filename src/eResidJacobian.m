
function JT = eResidJacobian(C,p,d)
JT = zeros(size(C,1),3);
for ii = 1:size(C,1);
    
    JT(ii,1) = -2*(p(3)-1)*(d(ii,1)-C(ii,1)+p(1)+(p(3)*C(ii,1))-(p(3)*p(1)));
    JT(ii,2) = -2*(p(3)-1)*(d(ii,2)-C(ii,2)+p(2)+(p(3)*C(ii,2))-(p(3)*p(2)));
    JT(ii,3) = 2*(C(ii,1)-p(1))*(d(ii,1)-C(ii,1)+p(1)+(p(3)*C(ii,1))-(p(3)*p(1))) + ...
        2*(C(ii,2)-p(2))*(d(ii,2)-C(ii,2)+p(2)+(p(3)*C(ii,2))-(p(3)*p(2)));
    
end