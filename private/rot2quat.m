function q = rot2quat(R)
    assert(all(size(R)==[3,3]))
    
    A =  [R(1,1)+R(2,2)+R(3,3),        -R(2,3)+R(3,2),       -R(3,1)+R(1,3),       -R(1,2)+R(2,1)
                -R(2,3)+R(3,2),  R(1,1)-R(2,2)-R(3,3),        R(2,1)+R(1,2),        R(3,1)+R(1,3)
                -R(3,1)+R(1,3),         R(2,1)+R(1,2), R(2,2)-R(1,1)-R(3,3),        R(3,2)+R(2,3)
                -R(1,2)+R(2,1),         R(3,1)+R(1,3),        R(3,2)+R(2,3), R(3,3)-R(1,1)-R(2,2)];
         
    [V,D] = eig(A - 3*eye(4));
    [~,ind] = min(abs(diag(D)));
    q = V(:,ind);
    
end