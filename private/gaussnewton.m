function [x, res, stop] = gaussnewton(fun, x0, jacobflag, maxitr, varargin)

    x      = x0;
    itr    = 0;
    stop   = 0;
    r      = [];
    tolx   = 1e-8;
    tolfun = 1e-8;
    while ~stop
        itr = itr + 1;
        
        if jacobflag
            [r, J] = feval(fun, x, varargin{:});
        else
            [J, r] = numericjacobian(fun, x, r, varargin{:});
        end
        
        x1 = x - J\r;
        r1 = feval(fun, x1, varargin{:});
        
        
        % stopping criteria
        dx = x1 - x;
        dr = r1 - r;
        if max(abs(dx)) < tolx
            stop = 1;
        elseif max(abs(dr)) < tolfun
            stop = 2;
        elseif max(abs(r1)) < tolfun
            stop = 3;
        elseif itr == maxitr
            stop = -1;
        end
        
        x = x1;
        r = r1;
        
        res = norm(r1)^2;
    end
    
end