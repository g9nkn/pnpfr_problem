# A Versatile Approach for Solving PnP, PnPf, and PnPfr Problems

&copy; 2020 NEC Corporation

This repository is an official MATLAB implementation of the paper "A Versatile Approach for Solving PnP, PnPf, and PnPfr Problems", ECCV2016 [(pdf)](https://jpn.nec.com/rd/people/docs/eccv2016_nakano.pdf) [(supp)](https://jpn.nec.com/rd/people/docs/eccv2016_appendix_nakano.pdf) [\[1\]](#reference).
A Gr&ouml;bner basis solver used in the code is generated by V. Larsson's automatic generator for polynomial solvers [\[2\]](#reference).  
For Japanese readers: 日本語によるPnPfr問題の解説は私の博士論文 [\[3\]](#reference) をご一読ください．

## License

This software is released under the NEC Corporation License.
See [LICENSE](https://github.com/g9nkn/pnpfr_problem/blob/main/LICENSE) before using the code. If you use this code, please cite the paper.

```bibtex
@inproceedings{nakano2016versatile,
  title={A versatile approach for solving PnP, PnPf, and PnPfr problems},
  author={Nakano, Gaku},
  booktitle={European Conference on Computer Vision},
  pages={338--352},
  year={2016},
  organization={Springer}
}
```

For commercial use, please contact Gaku Nakano \<g-nakano@nec.com\>.

## Usage

### Install

Copy `pnpfr_problem` folder and set a path to the folder by `addpath('pnpfr_problem')`. A subfolder `pnp_problem/private` is automatically linked.

### Function API

```
[R, t, err, f, k] = vpnpfr_nakano_eccv2016(pts3d, pts2d, problem, polishing)

INPUTS:
  pts3d - 3xN matrix of 3D points corresponding to the 2D points. (N >= 5)
          [x1, x2, ..., xN
           y1, y2, ..., yN
           z1, z2, ..., zN];
  pts2d - 2xN matrix of 2D points corresponding to the 3D points. (N >= 5)
          [u1, u2, ..., uN
           v1, v2, ..., vN];
          For PnP problem, each points are normalized by using the
          intrinsic paramters. For PnPf and PnPfr problems, each points
          are shifted by the image center.
  problem - which problem to be solved:
           'pnp'  : Perspective-n-point (PnP) problem 
           'pnpf' : PnP with unknown focal length
           'pnpfr': PnP with unknown focal length and unknown radial
                    distortion
  polishing - (optional) an integer to set the number of iterations of
              Gauss-Newton method for root polishing.
              If <= 0, the root polishing is not performed. (default: 0)

OUTPUS:
  R - 3x3xM rotation matrix. R(:,:,i) corresponds to t(:,i) and err(i).
  t - 3xM translation vector.
  err - 1xM algebraic cost of the solution.
  f - 1xM focal length. f(i) corresponds to R(:,:,i) and t(:,i).
      f=ones(1,M) for PnP problem.
  k - 3xM radial distortion. k(:,i) corresponds to R(:,:,i) and t(:,i).
      k=zeros(3,M) for PnP and PnPf problems.
```

As in Eqs. (1)-(3) in the paper, Fitzgibbon's division model is used to formulate the radial lens distortion. Please note that this is **NOT** compatible with [the OpenCV's definition](https://docs.opencv.org/master/d9/d0c/group__calib3d.html). Run `demo_vpnpfr_nakano_eccv2016.m` to understand how it works.

### Note

To reporoduce the Gr&ouml;bner basis solver for solving the first sub-problem, install [the Larsson's automatic generator](<http://people.inf.ethz.ch/vlarsson/misc/autogen_v0_5.zip>) and use the `pnpfr_problem/problem_vpnp.m` as follows:

```matlab
opts = default_options();
solv = generate_solver('vpnp', @problem_vpnp, opts);
```

As a result, a solver file `solver_vpnp.m` will be generated in `automatic_generator/solvers/`.

## Reference

1. Gaku Nakano, "A Versatile Approach for Solving PnP, PnPf, and PnPfr Problems," ECCV2016.  
[(pdf)](https://jpn.nec.com/rd/people/docs/eccv2016_nakano.pdf)
[(supp)](https://jpn.nec.com/rd/people/docs/eccv2016_appendix_nakano.pdf)

2. Viktor Larsson et al., "Efficient Solvers for Minimal Problems by Syzygy-based Reduction," CVPR 2017.  
Code: <http://people.inf.ethz.ch/vlarsson/misc/autogen_v0_5.zip>

3. 中野学，Perspective-n-Point問題とその派生問題に対する安定かつ高速な解法に関する研究，博士論文，筑波大学，2021年3月．
<https://jpn.nec.com/rd/people/docs/doctoral_thesis_nakano.pdf>

## Contributors

- Gaku Nakano, Central Research Laboratories, NEC Corporation.  
<g-nakano@nec.com>  
(ENG) <https://www.nec.com/en/global/rd/people/gaku_nakano.html>  
(JPN) <https://jpn.nec.com/rd/people/gaku_nakano.html>
