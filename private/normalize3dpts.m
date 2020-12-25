% shift 3D data so that variance = sqrt(2), mean = 0
function [Xnew, mean3d, var3d] = normalize3dpts( X )
	
	n = length(X);
	
    mean3d = mean(X, 2);
    Xnew = X - repmat(mean3d, 1, n);
    
    % variance (isotropic)
    var3d = sum( sqrt(sum( Xnew.^2 ) ) )/n;
	inv_var3d = 1/var3d;
	
%    % variance (anisotropic)
%    var3d = sqrt(sum(Xnew.^2, 2) )/n;
%    
%    var3d = diag(var3d);
%    inv_var3d = inv(var3d);
%    idx = find( inv_var3d==inf );
%    inv_var3d( idx ) = 1;
%    var3d(idx) = 1;
    
    Xnew = inv_var3d*Xnew;
	
return