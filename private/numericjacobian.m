% Approximate Jacobian by forward differences
function [J, f] = numericjacobian(fun, x, f, varargin)

    if isempty(f)
        f = feval(fun, x, varargin{:});
    end

    
    d  = 1e-7; 
    J  = zeros(length(f), length(x));
    xx = repmat(x, 1, length(x)) + d*diag( max(abs(x),d) );
    for j = 1:length(x)
        fp = feval(fun, xx(:,j), varargin{:});
        J(:,j) = (fp - f) / (xx(j,j) - x(j));
    end
    

end