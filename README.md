# Rational minimax approximation of matrix-valued functions

This folder contains the codes used to generate the numerical results in

L.-H. Zhang, Y.-N. Zhang, C. Zhang and S. Han,
*Rational minimax approximation of matrix-valued functions*, 2025.
Preprint: <https://arxiv.org/pdf/2508.06378>

A PDF of the manuscript is included as
[`mdlawson.pdf`](mdlawson.pdf).

The method introduced in that paper is implemented in `m_d_lawson.m`.
The scripts below reproduce the tables and figures:

| Result in the paper | Script |
| --- | --- |
| Table 8.1 | `run_table1_10times.m` |
| Table 8.2 | `run_table2_10times.m` |
| Figure 8.1 | `run_fig2_exam2_rat.m` |
| Figure 8.2 | `run_fig3_exam2_compare.m` |
| Figure 8.3 | `run_fig1_exam2_polys.m` |
| Figure 8.4 | `run_fig4_exam2_converg.m` |
| Figure 8.5 | `run_fig5_changheng.m` |

Before running these scripts, add `block_aaa-master`, `rktoolbox`, and
`matrix_fitting_toolbox_1` to the MATLAB path.

## Third-party software

block-AAA, RKFIT, and the Matrix Fitting Toolbox are **not** part of this
work. They are redistributed here only so that the numerical comparisons in
the paper can be reproduced. Copyright and license terms remain with their
original authors. Please cite the original papers and follow the original
licenses. Nothing in this repository transfers those rights.

### block-AAA

- Folder: `block_aaa-master/`
- Authors: Ion Victor Gosea and Stefan Güttel
- Paper: I. V. Gosea and S. Güttel, *Algorithms for the rational approximation of matrix-valued functions*, SIAM Journal on Scientific Computing, 2021. <https://arxiv.org/abs/2003.06410>
- Source code: <https://github.com/nla-group/block_aaa>
- License: MIT License, Copyright (c) 2019 nla-group. The full text is in `block_aaa-master/block_aaa-master/LICENSE`.

### RKFIT (Rational Krylov Toolbox)

- Folder: `rktoolbox/`
- RKFIT is part of the Rational Krylov Toolbox (RKToolbox), version 2.9 (14 July 2020), by the Rational Krylov Team at The University of Manchester.
- Copyright 2020 by the Rational Krylov Team, The University of Manchester. See `rktoolbox/Licence.txt` and `rktoolbox/Contents.m`.
- Website: <http://rktoolbox.org>
- The licence permits redistribution only if the original copyright notice and disclaimer are retained. The name of The University of Manchester and the names of its contributors may not be used to endorse products derived from this software without prior written permission.
- Reference: M. Berljafa and S. Güttel, *The RKFIT algorithm for nonlinear rational approximation*, SIAM Journal on Scientific Computing 39(5), A2049–A2071, 2017.

### Matrix Fitting Toolbox

- Folder: `matrix_fitting_toolbox_1/`
- Author: Bjørn Gustavsen, SINTEF Energy Research, Trondheim, Norway.
- The routine used in the experiments is `vectfit3.m` (Fast Relaxed Vector Fitting, version 1.0, 8 August 2008). Its header states that the software is **restricted to non-commercial use**.
- Website: <https://www.sintef.no/projectweb/vectorfitting/>
- Reference: B. Gustavsen and A. Semlyen, *Rational approximation of frequency domain responses by vector fitting*, IEEE Transactions on Power Delivery 14(3), 1052–1061, 1999.
- The PDF articles stored inside `matrix_fitting_toolbox_1/` belong to the original toolbox distribution. Their copyright remains with the original authors and publishers.
