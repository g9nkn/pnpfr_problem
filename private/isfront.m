function front = isfront(R, t, X, s)
    
    projective_depth = R(3,:)*X + t(3);
    front = nnz(projective_depth > 0)/size(X,2) > s; 
    
end