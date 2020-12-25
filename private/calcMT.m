function [M, T, A1, B1, A2, B2, D] = calcMT(X, pt2d)
	
    n = size(X, 2);

 	u = transpose( pt2d(1,:)./pt2d(3,:) );
 	v = transpose( pt2d(2,:)./pt2d(3,:) );
    X = X';
    
    uX = u.*X;
    vX = v.*X;
    
    A = [-vX, uX]/n;
    B = [-v, u]/n;
    
    T = -B\A;
    M = (A + B*T)' * (A + B*T);
    
    
    O  = zeros(n, 1);
    O3 = zeros(n, 3);
    l  =  ones(length(u), 1);
    A1 = [O3, -X,  vX;
           X, O3, -uX];
    B1 = [O, -l,  v;
          l,  O, -u];
    
    
    if nargout > 4
        r2 = u.^2 + v.^2;
        r4 = r2.^2;
        r6 = r2.*r4;
        D  = [r2, r4, r6];
        
        B2 = [O3, -D;
               D,  O3];
           
           
        O9 = zeros(n, 9);
        E  = [X(:,1).*D, X(:,2).*D, X(:,3).*D];
%         E  = [bsxfun(@times, X(:,1), D), bsxfun(@times, X(:,2), D), bsxfun(@times, X(:,3), D)];
        
        A2 = [O9,  -E;
               E,  O9];
    end
    
return