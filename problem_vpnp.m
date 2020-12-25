function [ eqs, data0, eqs_data ] = problem_vpnp(data0)
    
    if nargin < 1 || isempty(data0)
        M = randi(30097,6,6);
        M = M+M';
        data0 = M(:);
    end

    n  = 5;
    M  = reshape(data0,6,6);
    xx = create_vars(n);
    
    eqs = buildEqs(xx, M);
    
       
    if nargout == 3
        vars = create_vars(n+36);
        xx   = vars(1:n);
        data = vars(n+1:end);
        M    = reshape(data,6,6);
        
        eqs_data = buildEqs(xx, M);
    end

       
end


function eqs = buildEqs(xx, M)
    
    r1 = [xx(1); xx(2); xx(3)];
    r2 = [xx(4); xx(5);     1]; %linear constraint, R(2,3)=1

    % build polynomial equations
    A    = M(1:3,1:3)*r1  + M(1:3,4:6)*r2;
    B    = M(1:3,4:6)'*r1 + M(4:6,4:6)*r2;
    S1   = [     0, -r1(3),  r1(2); 
             r1(3),      0, -r1(1); 
            -r1(2),  r1(1),      0]; % [r1]_x
    S2   = [     0, -r2(3),  r2(2); 
             r2(3),      0, -r2(1); 
            -r2(2),  r2(1),      0]; % [r2]_x
    eq1  = S1*A + S2*B;              % Eq.(23)
    eq2  = r2'*A - r1'*B;            % Eq.(26)
    ceq1 = r1'*r1 - r2'*r2;          % Eq.(19)
    ceq2 = r1'*r2;                   % Eq.(20)
    eqs  = [eq1(:) 
            eq2
            ceq1
            ceq2];

end