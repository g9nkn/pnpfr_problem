function err = reprojerr(X, m, R, t, f, k)

    narginchk(4, 6);
    if nargin < 6, k = zeros(3,1); end
    if nargin < 5, f = 1; end
    
    
    Xc      = bsxfun(@plus, R*X, t);
    Xc(3,:) = 1/f * Xc(3,:);
    
    if size(m,1) == 3
        m = bsxfun(@times, m(1:2,:), 1./m(3,:));
    end
    r2 = sum(m(1:2,:).^2);
    r4 = r2.^2;
    r6 = r2.*r4;
    m  = bsxfun(@times, m(1:2,:), 1./(1 + k(1)*r2+k(2)*r4+k(3)*r6) );
    
    err = m(1:2,:) - bsxfun(@times, Xc(1:2,:), 1./Xc(3,:) );
    err = (norm(err, 'fro')^2) / (2*size(m,2));

end