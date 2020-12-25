function R = quat2mat(q)
% QUAT2MAT - Convert a quaternion to a 3x3 rotation matrix
%   

a = q(1); b = q(2); c = q(3); d = q(4);
R = zeros(3,3);
if isa(q, 'sym')
    R = sym(R);
end
R(1) = a^2+b^2-c^2-d^2;
R(2) = 2*(b*c+a*d);
R(3) = 2*(b*d-a*c);  
R(4) = 2*(b*c-a*d);
R(5) = a^2+c^2-b^2-d^2;
R(6) = 2*(c*d+a*b);
R(7) = 2*(b*d+a*c); 
R(8) = 2*(c*d-a*b);
R(9) = a^2+d^2-b^2-c^2;
end