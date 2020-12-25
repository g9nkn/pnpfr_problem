npts  = 100;
planar_scene = false;


% intrinsic parameters
w = 640;         % image size
h = 480;
cx = w/2;        % optical center (assumed to be image center)
cy = h/2;
cxy = [cx,cy]';
f = 400;         % focal length
K = [f, 0, cx
     0, f, cy
     0, 0,  1];
k = [-0.2/f^2    % radial distortion coeffs
            0
            0]; 

% distorted (oberved) points 
% generate points within [-w/2,-h/2] to [+w/2, +h/2]
u_d = w*rand(1, npts) - w/2;
v_d = h*rand(1, npts) - h/2;
pts2d_d = [u_d + cx
           v_d + cy
           ones(1,npts)];

% undistorted points
r2  = u_d.^2 + v_d.^2;
r4  = r2.^2;
r6  = r2.^3;
radial = ones(1, npts) + k(1)*r2 + k(2)*r4 + k(3)*r6;
u_u = u_d ./ radial;
v_u = v_d ./ radial;
pts2d_u = [u_u + cx
           v_u + cy
           ones(1,npts)];

       
% projective depth
min_depth = 5;
max_depth = 10;
if planar_scene
    % 3D plane: normal'*X+p = 0 where X=d*K\m
    theta  = (2*pi/4)*rand(1) - pi/4;
    phi    = pi/3*rand(1);
    normal = -[sin(theta)*cos(phi)
               sin(theta)*sin(phi) 
               cos(phi)];
    p = (max_depth - min_depth)*rand(1);
    d = - p ./ (normal'*(K\pts2d_u));
else
    d = min_depth + (max_depth - min_depth)*rand(1, npts);
end

       
% 3D points d*m = K*(R*X+t) -> X = R'*(d*K\m - t)
q = rand(4,1);
R = quat2rot( q/norm(q) );
t = min_depth/2 * rand(3,1);
pts3d = R' * (d.*(K\pts2d_u) - t);


% add image noise 
sigma = 1;
noisy_pts2d_u = pts2d_u(1:2,:) + sigma*randn(2,npts);
noisy_pts2d_d = pts2d_d(1:2,:) + sigma*randn(2,npts);


polishing = 0; % >=0, the number of iterations for root polishing

% solve pnp problem
pts2d = (noisy_pts2d_u - cxy) / f;
[R1, t1, err1] = vpnpfr_nakano_eccv2016(pts3d, pts2d, 'pnp', polishing);
disp('PnP Estimation errors:')
disp(['    R : ' num2str(min(calc_R_err(R, R1)))])
disp(['    t : ' num2str(min(calc_t_err(t, t1)))])

% solve pnpf problem
pts2d = noisy_pts2d_u - cxy;
[R2, t2, err2, f2] = vpnpfr_nakano_eccv2016(pts3d, pts2d, 'pnpf', polishing);
disp('PnPf Estimation errors:')
disp(['    R : ' num2str(min(calc_R_err(R, R2)))])
disp(['    t : ' num2str(min(calc_t_err(t, t2)))])
disp(['    f : ' num2str(min(calc_f_err(f, f2)))])

% solve pnpfr problem
pts2d = noisy_pts2d_d - cxy;
[R3, t3, err3, f3, k3] = vpnpfr_nakano_eccv2016(pts3d, pts2d, 'pnpfr', polishing);
disp('PnPfr Estimation errors:')
disp(['    R : ' num2str(min(calc_R_err(R, R3)))])
disp(['    t : ' num2str(min(calc_t_err(t, t3)))])
disp(['    f : ' num2str(min(calc_f_err(f, f3)))])
disp(['    k1: ' num2str(min(calc_k1_err(k(1,:), k3(1,:))))])


function R_err = calc_R_err(R_gt, R_est)
    for i=1:size(R_est,3)
        R_err(i) = norm(R_gt'*R_est(:,:,i) - eye(3), 'fro');
    end
end

function t_err = calc_t_err(t_gt, t_est)
    for i=1:size(t_est,2)
        t_err(i) = norm(t_gt - t_est(:,i)) / norm(t_gt);
    end
end

function f_err = calc_f_err(f_gt, f_est)
    for i=1:length(f_est)
        f_err(i) = abs(f_gt - f_est(i)) / f_gt;
    end
end

function k1_err = calc_k1_err(k1_gt, k1_est)
    for i=1:length(k1_est)
        k1_err(i) = abs(k1_gt - k1_est(i)) / abs(k1_gt);
    end
end
