function [R, t, f, err] = solver_vpnpf(X, m, polishing)
	
    narginchk(2, 3);
	
	%% calc coefficient matrices
	[M1, ~, A1, B1] = calcMT(X, m);	
    
	%% solve the 1st sub-problem by Groebner basis method
    R_tmp = solver_common_gb( M1 );
    
    %% solve the 2nd sub-problem by linear least-squares
    nsols = size(R_tmp,3);
    R     = zeros(3,3,nsols);
    t     = zeros(3,nsols);
    w     = zeros(1,nsols); % = 1/f, reciprocal of focal length
    err   = zeros(1,nsols);
    front = false(1,nsols);
    th    = 0.9; % be sure that 90% of points are in front of the camera

    T2 = -B1\A1;
    M2 = (A1+B1*T2)' * (A1+B1*T2);
    for i = 1:nsols
        R(:,:,i)       = R_tmp(:,:,i);
        [t(:,i), w(i)] = get_wt(R(:,:,i), M2, T2);
        err(i)         = algbraic_err(R(:,:,i), w(i), M2);
        front(i)       = isfront(R(:,:,i), t(:,i), X, th);        
    end

    % remove solutions that have "minus" focal length or 3d points don't
    % exist in front of the camera.
    valid = (w > 0) & (front > 0);
    R     = R(:,:,valid);
    t     = t(:,valid);
    w     = w(valid);
    err   = err(valid);

    %% perform root polishing
    if polishing > 0
        for i = 1:size(R,3)
            [R(:,:,i), w(i)] = rootpolishing(R(:,:,i), w(i), polishing, M2);
            r        = reshape(R(:,:,i)',9,1);
            W        = diag([1,1,1,1,1,1,w(i),w(i),w(i)]);
            invK     = diag([1,1,1/w(i)]);
            t(:,i)   = invK*T2*W*r;
            err(i)   = algbraic_err(R(:,:,i), w(i), M2);
            
            % flip if w is conveged to w < 0
            % diag([1,1,-w]*(R*X+t) = (diag([1,1,-w]*Rz) * Rz*(R*X+t) 
            if w(i) < 0
                Rz       = diag([-1,-1,1]);
                w(i)     = -w(i);
                R(:,:,i) = Rz*R(:,:,i);
                t(:,i)   = Rz*t(:,i);
            end
        end
    end

    f = 1./w;
	
end


function [t, w] = get_wt(R, M, T)
    r     = reshape(R', 9, 1);
    w     = - ( r(1:6)'*M(1:6,7:9)*r(7:9) ) / ( r(7:9)'*M(7:9,7:9)*r(7:9) );
    
    f    = 1/w;
    invK = diag([1, 1, f]);
    W    = diag([1,1,1,1,1,1,w,w,w]);

    t = invK * T * W * r;
end


function err = algbraic_err(R, w, M)
    W = diag([1,1,1,1,1,1,w,w,w]);
    r = reshape(R', 9, 1);
    
    err = norm( M*W*r )^2;
end

function [R, w, err] = rootpolishing(R, w, maxitr, M)
    
    q = rot2quat(R);
    x = [q; w];
    use_jacobian = false;
    [x, err] = gaussnewton(@objfunc, x, use_jacobian, maxitr, M);
    q = x(1:4);
    w = x(5);
    R = quat2rot(q/norm(q));
end


function [eqs, J] = objfunc(x, M)
    
    q = x(1:4);
    w = x(5);

    R    = quat2rot(q);    
    Rt   = R';
    r    = Rt(:);
    W    = diag([1,1,1,1,1,1,w,w,w]);
    dW   = diag([0,0,0,0,0,0,1,1,1]);
    WMWr = reshape(W*M*W*r,3,3)';

    A = Rt*WMWr;
    B = WMWr*Rt;

    eqs = [A([2,3,6]) - A([4,7,8]), ... 
           B([2,3,6]) - B([4,7,8]), ... 
           norm(q)^2 - 1, ...
           r'*(dW*M*W + W*M*dW)*r]';
       
    if nargout == 2
        J = jacobfunc(x, M);
    end
end


function J = jacobfunc(x, M)

    q = x(1:4);
    w = x(5);

    R    = quat2rot(q);
    Rt   = R';
    r    = Rt(:);
    W    = diag([1,1,1,1,1,1,w,w,w]);
    dW   = diag([0,0,0,0,0,0,1,1,1]);
    WMW  = W*M*W;
    WMWr = reshape(WMW*r,3,3)';
    G    = dW*M*W + W*M*dW;

    a = q(1); b = q(2); c = q(3); d = q(4);
    drdq = [  a,  b, -c, -d
             -d,  c,  b, -a
              c,  d,  a,  b
              d,  c,  b,  a
              a, -b,  c, -d
             -b, -a,  d,  c
             -c,  d, -a,  b
              b,  a,  d,  c
              a, -b, -c,  d];
   
   J = zeros(8, 5);
   for i = 1:4
       dR    = reshape(drdq(:,i),3,3)';
       dWMWr = reshape(WMW*drdq(:,i),3,3)';

       dA = dR'*WMWr + Rt*dWMWr;
       dB = WMWr*dR' + dWMWr*Rt;
       
       J(1:3, i) = dA([2,3,6]) - dA([4,7,8]);
       J(4:6, i) = dB([2,3,6]) - dB([4,7,8]);
       J(  7, i) = q(i);
       J(  8, i) = r'*(G+G')*drdq(:,i); 
   end
   J(:,1:4) = 2*J(:,1:4);
   
   dG    = reshape( G*r, 3, 3)';
   dAdw  = Rt*dG;
   dBdw  = dG*Rt;
   J(1:3, 5) = dAdw([2,3,6]) - dAdw([4,7,8]);
   J(4:6, 5) = dBdw([2,3,6]) - dBdw([4,7,8]);
   J(  7, 5) = 0;
   J(  8, 5) = 2*r'*dW*M*dW*r; 
   
end