function [Rfun, relgap, dk, err, pol, zers, bestW, info] = ...
    m_d_lawson(x1, bigF, nd, deg_vec, maxit)

% m-d-lawson is a method for solving rational approximation for matrix-valued
% functions. It is based on its dual form of this minimax approximation and
% the Lawson iteration is used to update the dual variable. 

% ---------------------------------------------------------
%
% Input:
%
%           x1                 sample points with length m
%           bigF           (m-by-n2) matrix [f_1,f_2,...,f_k]
%           nd               the degree of  denominator polynomial
%           deg_vec        degrees (n1,n2,...,nk) of numerators  p1,...,pk
%           maxit          maximun lawson's iteration step
% Output:
%
%          Rfun             cellarray with Rfun{1} = r1, ... , Rfun{ell} = r_ell
%                               each r_k = p_k/q ~ f_k is
%                           function handle with  degree (nk,nd)
%          pol              zeros of q
%          zers             cellarray with zers{k} is zeros of pk
%          dk               the sequence of the dual function values d2w
%          err              the record of maximun F-norm errors, i.e.,
%                               err(k) = MAX_j ||bigF(j,:)-bigR(j,:)||_{fro};
%          relgap           relative daulity gap at the computed approximant;
%                           generally, relgap < 0.01 means the computed 
%                           can be viewed as the minimax approximant;
%
%
% remark: this code is also suitable for the scalar-value function, as well
%         as the degree-specified duplexer problem
% ------------------------------------------------------------

% Reference:
% [1] L.-H. Zhang, Y.-N. Zhang, C. Zhang and S. Han, Rational minimax 
%     approximation of matrix-valued functions, 
%     2025, URL https://arxiv.org/pdf/2508.06378v2.

% [2] L.-H. Zhang and C. Zhang, Rational minimax approximations for 
%     matrix-valued functions:Existence, optimality and algorithms, 
%     URL https://arxiv.org/pdf/2607.22576.

% ----------------------------------------------------------------


    x1 = x1(:);
    [m, nResponses] = size(bigF);
    deg_vec = deg_vec(:).';

    validateattributes(x1, {'numeric'}, {'vector', 'finite', 'nonempty'}, ...
        mfilename, 'x1', 1);
    validateattributes(bigF, {'numeric'}, {'2d', 'finite', 'nonempty'}, ...
        mfilename, 'bigF', 2);
    validateattributes(nd, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'nd', 3);
    validateattributes(deg_vec, {'numeric'}, ...
        {'vector', 'integer', 'nonnegative', 'finite'}, ...
        mfilename, 'deg_vec', 4);
    validateattributes(maxit, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative', 'finite'}, ...
        mfilename, 'maxit', 5);

    if length(x1) ~= m
        error('m_d_lawson:SizeMismatch', ...
            'The number of sample points in x1 must match size(bigF,1).');
    end
    if length(deg_vec) ~= nResponses
        error('m_d_lawson:DegreeCount', ...
            'deg_vec must contain one numerator degree per response column.');
    end

    maxDegree = max([nd, deg_vec]);
    if maxDegree + 1 > m
        error('m_d_lawson:InsufficientSamples', ...
            'The requested maximum degree requires at least %d samples.', ...
            maxDegree + 1);
    end

    tol = 1e-2;
    numEvaluations = maxit + 1;
    w = ones(m, 1) / m;
    dk = nan(1, numEvaluations);
    err = nan(1, numEvaluations);
    gapHistory = nan(1, numEvaluations);

    bestEta = inf;
    bestK = 0;
    bestA = {};
    bestB = [];
    bestW = [];
    bestH = [];
    bestDual = NaN;
    stopReason = 'iteration limit';

    for kk = 1:numEvaluations
        [Q, H] = local_arnoldi_basis(x1, sqrt(w), maxDegree + 1);
        Qd = Q(:, 1:nd + 1);

        PFQ = zeros(m * nResponses, nd + 1, 'like', bigF);
        for j = 1:nResponses
            Qj = Q(:, 1:deg_vec(j) + 1);
            fQ = bigF(:, j) .* Qd;
            rows = (j - 1) * m + (1:m);
            PFQ(rows, :) = fQ - Qj * (Qj' * fQ);
        end

        [~, sigma, coeff] = svd(PFQ, 0);
        dualValue = sigma(end, end);
        b = coeff(:, end);
        q = Qd * b;

        Ak = cell(1, nResponses);
        bigR = zeros(size(bigF), 'like', bigF);
        for j = 1:nResponses
            Qj = Q(:, 1:deg_vec(j) + 1);
            aj = Qj \ (bigF(:, j) .* q);
            Ak{j} = aj(:);
            bigR(:, j) = (Qj * aj) ./ q;
        end

        perFrequencyError = sqrt(sum(abs(bigF - bigR).^2, 2));
        eta = max(perFrequencyError);
        if ~isfinite(eta)
            error('m_d_lawson:NonfiniteModel', ...
                'A nonfinite approximation was produced at iteration %d.', kk);
        end

        dk(kk) = dualValue;
        err(kk) = eta;
        gapHistory(kk) = local_relative_gap(eta, dualValue);

        if eta < bestEta
            bestEta = eta;
            bestA = Ak;
            bestB = b;
            bestW = w;
            bestH = H;
            bestDual = dualValue;
            bestK = kk;
        end

        if gapHistory(kk) < tol
            stopReason = 'relative-gap tolerance';
            break;
        end
        if eta < 2e-15
            stopReason = 'near-zero approximation error';
            break;
        end
        if kk < numEvaluations
            newW = w .* perFrequencyError;
            weightSum = sum(newW);
            if ~(isfinite(weightSum) && weightSum > 0)
                stopReason = 'zero Lawson weight update';
                break;
            end
            w = newW / weightSum;
        end
    end

    evaluatedIterations = kk;
    dk = dk(1:evaluatedIterations);
    err = err(1:evaluatedIterations);
    gapHistory = gapHistory(1:evaluatedIterations);

    if bestK == 0
        error('m_d_lawson:NoFiniteIterate', ...
            'No finite d-Lawson model was produced.');
    end

    Rfun = cell(1, nResponses);
    for j = 1:nResponses
        Rfun{j} = @(zz) zz_ratevalA(zz, bestA{j}, bestB, bestH);
    end

    if nd > 0
        [pol, ~] = zz_roots(bestH, bestB);
    else
        pol = [];
    end

    zers = cell(1, nResponses);
    for j = 1:nResponses
        if deg_vec(j) > 0
            [zers{j}, ~] = zz_roots(bestH, bestA{j});
        else
            zers{j} = [];
        end
    end

    selectedR = zeros(size(bigF), 'like', bigF);
    for j = 1:nResponses
        selectedR(:, j) = Rfun{j}(x1);
    end
    selectedPerFrequencyError = sqrt(sum(abs(bigF - selectedR).^2, 2));
    selectedEta = max(selectedPerFrequencyError);
    relgap = local_relative_gap(selectedEta, bestDual);

    info = struct();
    info.bestIteration = bestK;
    info.evaluatedIterations = evaluatedIterations;
    info.requestedLawsonUpdates = maxit;
    info.stopReason = stopReason;
    info.selectedError = selectedEta;
    info.selectedDualValue = bestDual;
    info.selectedRelgap = relgap;
    info.errorHistory = err;
    info.dualHistory = dk;
    info.gapHistory = gapHistory;
    info.perFrequencyError = selectedPerFrequencyError;
    % Retain the selected Arnoldi-basis representation for exact export to
    % pole-residue/state-space formats after the Lawson solve.
    info.arnoldiH = bestH;
    info.denominatorCoefficients = bestB;
    info.numeratorCoefficients = bestA;
    info.denominatorDegree = nd;
    info.numeratorDegrees = deg_vec;
    info.poles = pol(:);
    info.doesEnforceStability = false;
    info.doesEnforcePassivity = false;
    info.doesEnforceReciprocity = false;
    info.isStableInCurrentVariable = isempty(pol) || all(real(pol) < 0);
end

function gap = local_relative_gap(eta, dualValue)
    if eta <= eps
        gap = 0;
    else
        gap = abs(eta - dualValue) / eta;
    end
end

function [Q, H] = local_arnoldi_basis(nodes, startVector, numColumns)
    m = length(nodes);
    Q = zeros(m, numColumns, 'like', nodes);
    H = zeros(numColumns, max(numColumns - 1, 0), 'like', nodes);

    startNorm = norm(startVector);
    if ~(isfinite(startNorm) && startNorm > 0)
        error('m_d_lawson:InvalidWeight', ...
            'Lawson weights do not define a valid initial Arnoldi vector.');
    end
    Q(:, 1) = startVector / startNorm;

    for k = 1:numColumns - 1
        v = nodes .* Q(:, k);
        for pass = 1:2
            for j = 1:k
                h = Q(:, j)' * v;
                H(j, k) = H(j, k) + h;
                v = v - Q(:, j) * h;
            end
        end

        hNext = norm(v);
        breakdownTol = 100 * eps(max(1, norm(nodes, inf)));
        if ~(isfinite(hNext) && hNext > breakdownTol)
            error('m_d_lawson:ArnoldiBreakdown', ...
                ['Arnoldi breakdown at basis column %d. The requested ', ...
                 'degree is incompatible with these samples and weights.'], k + 1);
        end
        H(k + 1, k) = hNext;
        Q(:, k + 1) = v / hNext;
    end
end
