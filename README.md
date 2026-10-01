# Rational minimax approximation of matrix-valued functions

MATLAB code for the numerical experiments in

> L.-H. Zhang, Y.-N. Zhang, C. Zhang, and S. Han,
> *Rational minimax approximation of matrix-valued functions*.

An earlier description of the **m-d-Lawson** method is the preprint
[arXiv:2508.06378](https://arxiv.org/pdf/2508.06378).
This repository reproduces the updated manuscript [`mdlawson.pdf`](mdlawson.pdf).
The codes, tables, and figures here correspond to `mdlawson.pdf`, not to that earlier arXiv version.

## Requirements

- MATLAB
- The three comparison packages shipped in this repository (see [Third-party software](#third-party-software))

Add them to the path before running the scripts:

```matlab
addpath(genpath('block_aaa-master'))
addpath('rktoolbox')
addpath('matrix_fitting_toolbox_1')
```

## Reproducing the paper

| Result in `mdlawson.pdf` | Script |
| --- | --- |
| Table 8.1 | `run_table1_10times` |
| Table 8.2 | `run_table2_10times` |
| Figure 8.1 | `run_fig2_exam2_rat` |
| Figure 8.2 | `run_fig3_exam2_compare` |
| Figure 8.3 | `run_fig1_exam2_polys` |
| Figure 8.4 | `run_fig4_exam2_converg` |
| Figure 8.5 | `run_fig5_changheng` |

`run_table1_10times` and `run_table2_10times` each average ten runs and write `table1_10times.tex` and `table2_10times.tex`. `run_table1` and `run_table2` return the errors and the runtime; they do not write a TeX file.

## How the code is called

`m_d_lawson` is the solver. The experiment scripts build the samples, call the solver, and then evaluate the error. Comparison methods are called beside it. They are not called from inside `m_d_lawson`.

```text
run_table1_10times
  └── run_table1
        ├── m_d_lawson
        ├── zz_v_aaa_lawson
        ├── block_aaa
        ├── baryfit6          (RKFIT)
        ├── loewfit3          (Loewner)
        └── baryfit5          (vector fitting)
              └── comput_errors_mrat / comput_errors_mrat2

run_table2_10times
  └── run_table2
        └── (same six calls as run_table1)

run_fig1_exam2_polys ── m_d_lawson
run_fig2_exam2_rat   ── m_d_lawson
run_fig3_exam2_compare
  └── same six calls as run_table1
run_fig4_exam2_converg
  ├── zz_v_aaa_lawson
  ├── m_d_lawson
  └── zz_v_dlawson      (earlier implementation, one comparison only)
run_fig5_changheng
  ├── zz_v_aaa_lawson
  └── m_d_lawson
```

Call the solver as

```matlab
Rfun = m_d_lawson(x1, bigF, nd, deg_vec, maxit);
```

| Argument | Meaning |
| --- | --- |
| `x1` | Sample nodes, length `m` |
| `bigF` | Sampled matrix entries, size `m` by the number of entries |
| `nd` | Degree of the common denominator |
| `deg_vec` | Numerator degree of each entry |
| `maxit` | Maximum number of Lawson updates |

`Rfun` is a cell of function handles, one rational function per matrix entry.

Inside each Lawson step, `m_d_lawson` does the following.

1. `local_arnoldi_basis` builds a Vandermonde–Arnoldi basis from the current weights.
2. `svd` computes the smallest singular value of the linearized dual residual and the denominator coefficients.
3. A least-squares solve produces the numerator coefficients and the maximum error.
4. `local_relative_gap` tests the duality gap. If the iteration continues, the weights are updated by `w ← w .* error` and then normalized.

After the loop, the selected coefficients are turned into a rational function:

- `zz_ratevalA` evaluates each entry.
- `zz_roots` computes the poles of the denominator and the zeros of each numerator.

The table scripts then call `comput_errors_mrat` or `comput_errors_mrat2`. Those error evaluations are not included in the reported runtime.

## Third-party software

block-AAA, RKFIT, and the Matrix Fitting Toolbox are not part of this work. They are included only so that the numerical comparisons can be reproduced. Copyright and license terms remain with their authors. Please cite the original papers and follow the original licenses.

### block-AAA

- Folder: `block_aaa-master/`
- Authors: Ion Victor Gosea and Stefan Güttel
- Paper: I. V. Gosea and S. Güttel, *Algorithms for the rational approximation of matrix-valued functions*, SIAM Journal on Scientific Computing, 2021. <https://arxiv.org/abs/2003.06410>
- Source: <https://github.com/nla-group/block_aaa>
- License: MIT License, Copyright (c) 2019 nla-group. See `block_aaa-master/block_aaa-master/LICENSE`.

### RKFIT (Rational Krylov Toolbox)

- Folder: `rktoolbox/`
- RKFIT is part of the Rational Krylov Toolbox (RKToolbox), version 2.9 (14 July 2020), by the Rational Krylov Team at The University of Manchester.
- Copyright 2020 by the Rational Krylov Team, The University of Manchester. See `rktoolbox/Licence.txt` and `rktoolbox/Contents.m`.
- Website: <http://rktoolbox.org>
- Redistribution must keep that copyright notice. The name of The University of Manchester and the names of its contributors may not be used to endorse derived products without prior written permission.
- Reference: M. Berljafa and S. Güttel, *The RKFIT algorithm for nonlinear rational approximation*, SIAM Journal on Scientific Computing 39(5), A2049–A2071, 2017.

### Matrix Fitting Toolbox

- Folder: `matrix_fitting_toolbox_1/`
- Author: Bjørn Gustavsen, SINTEF Energy Research, Trondheim, Norway.
- The routine used here is `vectfit3.m` (Fast Relaxed Vector Fitting, version 1.0, 8 August 2008). Its header states that the software is **restricted to non-commercial use**.
- Website: <https://www.sintef.no/projectweb/vectorfitting/>
- Reference: B. Gustavsen and A. Semlyen, *Rational approximation of frequency domain responses by vector fitting*, IEEE Transactions on Power Delivery 14(3), 1052–1061, 1999.
- The PDF articles in `matrix_fitting_toolbox_1/` belong to the original toolbox. Their copyright remains with the original authors and publishers.

## Citation

If you use the m-d-Lawson code, cite the updated manuscript `mdlawson.pdf`. The preprint [arXiv:2508.06378](https://arxiv.org/pdf/2508.06378) describes an earlier version of the method.
