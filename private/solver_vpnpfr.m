function [R, t, f, k, err] = solver_vpnpfr(X, m, polishing)
	
    narginchk(2, 3);
	
	%% calc coefficient matrices
	[M, T, A1, B1, A2, B2, d] = calcMT(X, m);	
    
	%% solve the 1st sub-problem by Groebner basis method
    R_tmp = solver_common_gb( M );

    %% solve the 2nd sub-problem by linear least-squares
    nsols = size(R_tmp,3);
    R     = zeros(3,3,nsols);
    t     = zeros(3,nsols);
    w     = zeros(1,nsols); % = 1/f, reciprocal of focal length
    k     = zeros(3,nsols);
    err   = zeros(1,nsols);
    front = false(1,nsols);
    th    = 0.9; % be sure that 90% of points are in front of the camera

    A1A1 = A1'*A1;
    A1A2 = A1'*A2;
    A1B1 = A1'*B1;
    A1B2 = A1'*B2;
    A2A2 = A2'*A2;
    A2B1 = A2'*B1;
    A2B2 = A2'*B2;
    B1B1 = B1'*B1;
    B1B2 = B1'*B2;
    B2B2 = B2'*B2;

    for i = 1:nsols
        R(:,:,i)               = R_tmp(:,:,i);
        [t(:,i), w(i), k(:,i)] = find_twk(R(:,:,i), T, X, m, d);
        err(i)                 = algbraic_err(R(:,:,i), t(:,i), w(i), k(:,i), A1A1, A1A2, A1B1, A1B2, A2A2, A2B1, A2B2, B1B1, B1B2, B2B2);
        front(i)               = isfront(R(:,:,i), t(:,i), X, th);
    end

    % remove solutions that have "minus" focal length or 3d points don't
    % exist in front of the camera.
    valid = (w > 0) & (front > 0);
    R     = R(:,:,valid);
    t     = t(:,valid);
    w     = w(valid);
    k     = k(:,valid);
    err   = err(valid);

    %% perform root polishing
    if polishing > 0
        for i = 1:size(R,3)
            [R(:,:,i), t(:,i), w(i), k(:,i)] = rootpolishing(R(:,:,i), w(i), k(:,i), polishing, A1A1, A1B1, A1A2, A1B2, B1B1, B1B2, A2A2, A2B1, A2B2, B2B2);
            err(i)  = algbraic_err(R(:,:,i), t(:,i), w(i), k(:,i), A1A1, A1A2, A1B1, A1B2, A2A2, A2B1, A2B2, B1B1, B1B2, B2B2);
            
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


function [t, w, k] = find_twk(R, T, X, m, d)
    r   = reshape(R', 9, 1);
    t12 = T * r(1:6);
    
    xc = ( R(1,:) * X + t12(1) )';
    yc = ( R(2,:) * X + t12(2) )';
    zc = ( R(3,:) * X )';
    u = m(1,:)';
    v = m(2,:)';
    
    
    % Eq.(29)
    L = [ v,  v.*zc, -yc.*d
         -u, -u.*zc,  xc.*d];
     
    g = [ yc
         -xc];
    
    t3wk = L\g;
    
    w = t3wk(2);
    k = t3wk(3:end);
    t = [t12; t3wk(1)/w];
    
end



function err = algbraic_err(R, t, w, k, A1A1, A1A2, A1B1, A1B2, A2A2, A2B1, A2B2, B1B1, B1B2, B2B2)
    K = diag([1,1,w]);
    W = diag([ones(1,6), w*ones(1,3)]);
    r = reshape(R', 9, 1);
    r = W*r;
    t = K*t;
    
%     % kron() is very slow
%     y = [kron(r(1:3),k)
%          kron(r(4:6),k)];
%     z = [t(1)*k
%          t(2)*k];

    % faster version
    y = zeros(18,1);
    y(1:3) = R(1,1)*k;
    y(4:6) = R(1,2)*k;
    y(7:9) = R(1,3)*k;
    y(10:12) = R(2,1)*k;
    y(13:15) = R(2,2)*k;
    y(16:18) = R(2,3)*k;

    z = zeros(6,1);
    z(1:3) = t(1)*k;
    z(4:6) = t(2)*k;

    err = r'*(A1A1*r + 2*A1A2*y + 2*A1B1*t + 2*A1B2*z) + ...
          y'*(A2A2*y + 2*A2B1*t + 2*A2B2*z) + ...
          t'*(B1B1*t + 2*B1B2*z) + ...
          z'*B2B2*z;
end



function [R, t, w, k, err] = rootpolishing(R, w, k, maxitr, varargin)

    q = rot2quat(R);
    x = [q; w; k];
    use_jacobian = false;
    [x, err, stop] = gaussnewton(@objfunc, x, use_jacobian, maxitr, varargin{:});
    q = x(1:4);
    w = x(5);
    k = x(6:8);
    
    
    
    % recover t
    R = quat2mat(q/norm(q));
    r = reshape(R',9,1);
    
    [A1, B1, A2, B2] = varargin{:};
    [C, D] = calcCD(k, w, A1, B1, A2, B2); % Eqs.(C-1) and (C-2)
    t = -(D\C)*r; % Eq.(31)
    
end


function err = objfunc(x, A1A1, A1B1, A1A2, A1B2, B1B1, B1B2, A2A2, A2B1, A2B2, B2B2)

    q = x(1:4);
    w = x(5);
    k = x(6:8);
    
    K = diag([1,1,w]);
    W = diag([1,1,1,1,1,1,w,w,w]);
    
    Y = zeros(18,9);
    Y(  1:3,1) = k;
    Y(  4:6,2) = k;
    Y(  7:9,3) = k;
    Y(10:12,4) = k;
    Y(13:15,5) = k;
    Y(16:18,6) = k;
    
    Z = zeros(6,3);
    Z(  1:3,1) = k;
    Z(  4:6,2) = k;
    
    CC = W*A1A1*W  + W*A1A2*Y  + (W*A1A2*Y)' + Y'*A2A2*Y;
    DD = K*B1B1*K  + K*B1B2*Z  + (K*B1B2*Z)' + Z'*B2B2*Z;
    DC = K*A1B1'*W + K*A2B1'*Y + Z'*A1B2'*W  + Z'*A2B2'*Y;
    

    R = quat2mat(q/norm(q));
    r = reshape(R',9,1);
    
    
    t = -(DD\DC)*r;
    
    E1 = zeros(9,2);
    E1(7:9,1) = r(7:9);
    E1(1:6,2) = r(1:6);
    
    E2 = zeros(3,2);
    E2(3,1)   = t(3);
    E2(1:2,2) = t(1:2);
    
    % gradient of objective function w.r.t w
    dedw = [1, 0]*( (E1'*A1A1*E1 + E1'*A1B1*E2 + (E1'*A1B1*E2)' + E2'*B1B1*E2)*[w;1] + (E1'*A1A2 + E2'*A2B1')*Y*r + (E1'*A1B2 + E2'*B1B2)*Z*t );
    
    I = eye(3);
    F1 = [R(1,1)*I
          R(1,2)*I
          R(1,3)*I
          R(2,1)*I
          R(2,2)*I
          R(2,3)*I];
    F2 = [t(1)*I
          t(2)*I];
    
    
    % gradient of objective function w.r.t k
    dedk = (F1'*A2A2*F1 + F1'*A2B2*F2 + (F1'*A2B2*F2)' + F2'*B2B2*F2)*k + (F1'*A1A2' + F2'*A1B2') *W*r + (F1'*A2B1 + F2'*B1B2') *K*t;
    
    
    M = CC - DC'*(DD\DC);
    Mr = reshape(M*r, 3, 3)';
    Rt = R';
    
    F1 = Rt*Mr;
    F2 = Mr*Rt;
    err = [F1([2,3,6]) - F1([4,7,8]), ... % Eq.(34)
           F2([2,3,6]) - F2([4,7,8]), ... % Eq.(35)
           norm(q)^2 - 1,...              % coincides with Eqs.(36) and (37)
           dedw,...                       % Eq.(38)
           dedk']';                       % Eq.(39)
end



% Eqs.(C-1) and (C-2)
function [C, D, A1W, B1K, A2Y, B2Z] = calcCD(k, w, A1, B1, A2, B2)
    K = diag([1,1,w]);
    W = diag([1,1,1,1,1,1,w,w,w]);
    
%   % kron() is very slow
%     Y = [kron(eye(6),k), zeros(18,3)];
%     Z = [kron(eye(2),k), zeros(6,1)];
    Y = zeros(18,9);
    Y(  1:3,1) = k;
    Y(  4:6,2) = k;
    Y(  7:9,3) = k;
    Y(10:12,4) = k;
    Y(13:15,5) = k;
    Y(16:18,6) = k;
    
    Z = zeros(6,3);
    Z(  1:3,1) = k;
    Z(  4:6,2) = k;
    

    A1W = A1*W;
    B1K = B1*K;
    A2Y = A2*Y;
    B2Z = B2*Z;
    
    C = A1W + A2Y;
    D = B1K + B2Z;
    
end