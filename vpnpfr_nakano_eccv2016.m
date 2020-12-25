%% A MATLAB code for solving PnP, PnPf, and PnPfr problems by G. Nakano [1].
% Copyright (c) 2020 NEC Corporation
% This software is released under the NEC Corporation License. See License.txt.
% For commercial use, please contact Gaku Nakano <g-nakano@nec.com>.
%
% USAGE:
%   [R, t, err, f, k] = vpnpfr_nakano_eccv2016(pts3d, pts2d, problem, polishing)
%
% INPUTS:
%   pts3d - 3xN matrix of 3D points corresponding to the 2D points. (N >= 5)
%           [x1, x2, ..., xN
%            y1, y2, ..., yN
%            z1, z2, ..., zN];
%   pts2d - 2xN matrix of 2D points corresponding to the 3D points. (N >= 5)
%           [u1, u2, ..., uN
%            v1, v2, ..., vN];
%           For PnP problem, each points are normalized by using the
%           intrinsic paramters. For PnPf and PnPfr problems, each points
%           are shifted by the image center.
%   problem - which problem to be solved:
%            'pnp'  : Perspective-n-point problem 
%            'pnpf' : PnP with unknown focal length
%            'pnpfr': PnP with unknown focal length and unknown radial
%                     distortion
%             The common solver was geenrated by using [2].
%   polishing - (optional) an integer to set the number of iterations of
%               Gauss-Newton method for root polishing.
%               If <= 0, the root polishing is not performed. (default: 0)
%
% OUTPUS:
%   R - 3x3xM rotation matrix. R(:,:,i) corresponds to t(:,i) and err(i).
%   t - 3xM translation vector.
%   err - 1xM algebraic cost of the solution.
%   f - 1xM focal length. f(i) corresponds to R(:,:,i) and t(:,i).
%       f=ones(1,M) for PnP problem.
%   k - 3xM radial distortion. k(:,i) corresponds to R(:,:,i) and t(:,i).
%       k=zeros(3,M) for PnP and PnPf problems.
%
% REFERENCE:
%   [1] Gaku Nakano, "A versatile approach for solving PnP, PnPf, and PnPfr
%       problems," ECCV2016.
%   [2] V. Larsson et al., "Efficient Solvers for Minimal Problems by
%       Syzygy-based Reduction," CVPR 2017.
%       http://people.inf.ethz.ch/vlarsson/misc/autogen_v0_5.zip
function [R, t, err, f, k] = vpnpfr_nakano_eccv2016(pts3d, pts2d, problem, polishing)

    %% check arguments
    narginchk(3,4);
    if nargin < 4
        polishing = 0;
    end
    assert( size(pts3d,1)==3 )
    assert( size(pts2d,1)==2 )
    assert( size(pts3d,2)==size(pts2d,2) && size(pts3d,2) >= 5 )
    
    %% random rotation to avoid degeneracy in case of r23=0
    % R*X = R * (R_rand'*R_rand) * X = (R*R_rand')*(R_rand*X)
    q      = rand(4,1);
    q      = q/norm(q);
    R_rand = quat2rot( q );
    pts3d  = R_rand * pts3d;
    
    %% normalize 3d points
    [pts3d_n, mean3d, scale3d] = normalize3dpts(pts3d);   
    
    %% normalize 2d points for pnpf or pnpfr
    if contains(problem,'pnpf')
        scale2d = max( abs(pts2d(:)) );
    else
        scale2d = 1;
    end
    pts2d_n      = pts2d / scale2d; 
    pts2d_n(3,:) = 1; % make 2d points as homogeneus coordinages

    %% solve the 1st and the 2nd sub-problems
    switch lower(problem)
        case 'pnp'
            [R, t, err] = solver_vpnp(pts3d_n, pts2d_n, polishing);
            f = ones(1,size(R,3))/scale2d;
            k = zeros(3,size(R,3));
        case 'pnpf'
            [R, t, f, err] = solver_vpnpf(pts3d_n, pts2d_n, polishing);
            k = zeros(3,size(R,3));
        case 'pnpfr'
            [R, t, f, k, err] = solver_vpnpfr(pts3d_n, pts2d_n, polishing);
        otherwise
            error('''problem'' is wrong. Check valid options by ''help vpnpfr_nakano_eccv2016''');
    end
    
    %% recover R and t in the original world coordinates
    for i=1:size(R,3)
        t(:,i)   = scale3d*t(:,i) - R(:,:,i)*mean3d;
        R(:,:,i) = R(:,:,i) * R_rand;
    end
    
    %% de-normalization
    f = f * scale2d;
    k = k ./ [scale2d^2
              scale2d^4
              scale2d^6];
end