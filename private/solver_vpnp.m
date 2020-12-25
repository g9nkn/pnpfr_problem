function [R, t, err] = solver_vpnp(X, m, polishing)
	
    narginchk(2, 3);
	
	%% calc coefficient matrices
	[M1, ~, A1, B1] = calcMT(X, m);	
    
	%% solve the 1st sub-problem by Groebner basis method
    R_tmp = solver_common_gb( M1 );
    

    %% solve the 2nd sub-problem by linear least-squares
    nsols = size(R_tmp,3);
    R     = zeros(3,3,nsols);
    t     = zeros(3,nsols);
    err   = zeros(1,nsols);
    
    front = false(1,nsols);
    th    = 0.9; % be sure that 90% of points are in front of the camera
    
    T2 = -B1\A1;
    M2 = (A1+B1*T2)' * (A1+B1*T2);
    for i = 1:nsols            
        R(:,:,i) = R_tmp(:,:,i);
        r        = reshape(R(:,:,i)',9,1);
        t(:,i)   = T2*r;
        err(i)   = r'*M2*r;        
        front(i) = isfront(R(:,:,i), t(:,i), X, th);
    end
    R   = R(:,:,front);
    t   = t(:,front);
    err = err(front);
    
    %% perform root polishing
    if polishing > 0
        for i=1:size(R,3)
            R(:,:,i) = rootpolishing(R(:,:,i), M2, polishing);
            r        = reshape(R(:,:,i)',9,1);
            t(:,i)   = T2*r;
            err(i)   = r'*M2*r; 
        end
    end
	

end


function [R, err] = rootpolishing(R, M, maxitr)
    
    q = rot2quat(R);
    use_jacobian = true;
    [q, err, ~] = gaussnewton(@objfunc, q, use_jacobian, maxitr, M);
    R = quat2rot(q);

end




function [eqs, J] = objfunc(q, M)
    
    R = quat2rot(q);
    Rt    = R';
    r     = Rt(:);
    matMr = reshape(M*r,3,3)';

    A = Rt*matMr;
    B = matMr*Rt;

    eqs = [A([2,3,6]) - A([4,7,8]), ...
           B([2,3,6]) - B([4,7,8]), ...
           norm(q)^2 - 1]';       
       
    if nargout == 2
        J = jacobfunc(q, M);
    end
end


function J = jacobfunc(q, M)
    a = q(1); b = q(2); c = q(3); d = q(4);
    R = quat2rot(q);
    
    r     = reshape(R', 9, 1);
    matMr = reshape(M*r,3,3)';
    
    
    drdq = [  a,  b, -c, -d
             -d,  c,  b, -a
              c,  d,  a,  b
              d,  c,  b,  a
              a, -b,  c, -d
             -b, -a,  d,  c
             -c,  d, -a,  b
              b,  a,  d,  c
              a, -b, -c,  d];
   
   J = zeros(7, 4);
   for i = 1:4
       dR  = reshape(drdq(:,i),3,3)';
       dMr = reshape(M*drdq(:,i),3,3)';

       dA = dR'*matMr + R'*dMr;
       dB = matMr*dR' + dMr*R';
       
       J(1:3, i) = dA([2,3,6]) - dA([4,7,8]);
       J(4:6, i) = dB([2,3,6]) - dB([4,7,8]);
       J(  7, i) = q(i);
       
   end
   J = 2*J;
   
end