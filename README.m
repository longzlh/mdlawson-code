% This folder contains all codes to generate the results for the paper:
%
% L.-H. Zhang, Y.-N. Zhang, C. Zhang and S. Han, Rational minimax
% approximation of matrix-valued functions, 2025,
% URL https://arxiv.org/pdf/2508.06378
%
% A PDF of the manuscript is included as
% mdlawson.pdf.
%
% NOTE:
%
% 1. The main algorithm is m_d_lawson.m.
%
% 2. To generate the results, first add these three third-party packages
%    to the MATLAB path:
%    block_aaa-master, rktoolbox, matrix_fitting_toolbox_1.
%    They are NOT part of this work. They are included only so that the
%    numerical comparisons can be reproduced. Copyright and licenses
%    remain with the original authors.
%
%    block-AAA
%      Folder: block_aaa-master
%      Authors: Ion Victor Gosea and Stefan Guettel
%      Paper: Algorithms for the rational approximation of matrix-valued
%        functions, SIAM J. Sci. Comput., 2021.
%        https://arxiv.org/abs/2003.06410
%      Source: https://github.com/nla-group/block_aaa
%      License: MIT, Copyright (c) 2019 nla-group.
%        See block_aaa-master/block_aaa-master/LICENSE.
%
%    RKFIT (Rational Krylov Toolbox)
%      Folder: rktoolbox
%      RKToolbox version 2.9 (14 July 2020), Rational Krylov Team,
%      The University of Manchester. http://rktoolbox.org
%      Copyright 2020. See rktoolbox/Licence.txt.
%      Redistribution must keep that copyright notice. The name of
%      The University of Manchester and of its contributors may not be
%      used to endorse derived products without prior written permission.
%      Reference: M. Berljafa and S. Guettel, The RKFIT algorithm for
%      nonlinear rational approximation, SIAM J. Sci. Comput. 39 (2017).
%
%    Matrix Fitting Toolbox
%      Folder: matrix_fitting_toolbox_1
%      Author: Bjorn Gustavsen, SINTEF Energy Research.
%      vectfit3.m is Fast Relaxed Vector Fitting, version 1.0
%      (8 August 2008), and is RESTRICTED to NON-COMMERCIAL use.
%      Website: https://www.sintef.no/projectweb/vectorfitting/
%      Reference: B. Gustavsen and A. Semlyen, Rational approximation of
%      frequency domain responses by vector fitting, IEEE Trans. Power
%      Delivery 14 (1999).
%      PDF files inside this folder belong to the original toolbox.
%      Their copyright remains with the original authors and publishers.
%
% 3. Correspondence to the paper:
%    Table 8.1  --- run_table1_10times.m
%    Table 8.2  --- run_table2_10times.m
%    Figure 8.1 --- run_fig2_exam2_rat.m
%    Figure 8.2 --- run_fig3_exam2_compare.m
%    Figure 8.3 --- run_fig1_exam2_polys.m
%    Figure 8.4 --- run_fig4_exam2_converg.m
%    Figure 8.5 --- run_fig5_changheng.m
