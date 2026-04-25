function results = solve_pipe_network(net)
%SOLVE_PIPE_NETWORK Hardy-Cross solver for closed-loop pipe networks.
%   This version matches the loop-correction workflow used in many
%   Fluid Mechanics courses (like the Daugherty/Ingersoll style example).
%
% Required fields
%   net.Q0        [nPipe x 1] initial flow guess in each branch
%   net.C         [nLoop x nPipe] loop-incidence matrix
%                 +1 branch aligned with loop direction
%                 -1 branch opposite to loop direction
%                  0 branch not in that loop
%
% Either provide constant K directly, or provide geometry for automatic K:
%   Option A (constant loss coefficients):
%       net.K     [nPipe x 1], where dh_i = K_i*Q_i*abs(Q_i)
%   Option B (Darcy-Weisbach updated each iteration):
%       net.L, net.D, net.eps, net.nu, net.g
%
% Optional
%   net.maxIter   default 200
%   net.tol       default 1e-8 (max |dQ_loop|)
%   net.nPrint    default 4 (rows captured for Table-3 style output)
%   net.pathList  cell array of candidate paths, each path is vector of
%                 branch indices used to compute total path head loss.
%
% Outputs
%   results.iterationTableQ  first nPrint iterations of Qi
%   results.Q                converged branch flows
%   results.dh               converged branch head loss per branch
%   results.Re               Reynolds number per branch (if D and nu given)
%   results.isTurbulent      Re >= 4000 (if Re is available)
%   results.highestLossPath  best path from pathList (if provided)

Q = net.Q0(:);
C = net.C;
[nLoop, nPipe] = size(C);

if nPipe ~= 22
    warning('Assignment expects 22 branches; C currently has %d branches.', nPipe);
end

if numel(Q) ~= nPipe
    error('Q0 size (%d) must match number of columns in C (%d).', numel(Q), nPipe);
end

maxIter = getFieldOr(net, 'maxIter', 200);
tol = getFieldOr(net, 'tol', 1e-8);
nPrint = getFieldOr(net, 'nPrint', 4);

Qhist = zeros(min(nPrint, maxIter), nPipe);
dqHist = zeros(maxIter, nLoop);

for k = 1:maxIter
    if nLoop == 0
        error('Loop-incidence matrix net.C is empty. Define loops from your network diagram.');
    end

    K = branchK(net, Q);

    for j = 1:nLoop
        idx = find(C(j,:) ~= 0);
        sgn = C(j, idx).';

        dhLoopTerms = sgn .* (K(idx) .* Q(idx) .* abs(Q(idx)));
        hLoop = sum(dhLoopTerms);

        denom = sum(2 .* K(idx) .* abs(Q(idx)));
        if denom == 0
            dq = 0;
        else
            dq = -hLoop / denom;
        end
        dqHist(k,j) = dq;

        Q(idx) = Q(idx) + sgn .* dq;
    end

    if k <= nPrint
        Qhist(k,:) = Q.';
    end

    if max(abs(dqHist(k,:))) < tol
        break;
    end
end

iterUsed = k;
K = branchK(net, Q);
dh = K .* Q .* abs(Q);

results.iterations = iterUsed;
results.iterationTableQ = Qhist(1:min(nPrint,iterUsed), :);
results.Q = Q;
results.K = K;
results.dh = dh;
results.dqPerLoop = dqHist(1:iterUsed, :);

if isfield(net, 'D') && isfield(net, 'nu')
    A = pi .* (net.D(:).^2) ./ 4;
    V = abs(Q) ./ A;
    Re = V .* net.D(:) ./ net.nu;
    results.Re = Re;
    results.isTurbulent = Re >= 4000;
else
    results.Re = [];
    results.isTurbulent = [];
end

if isfield(net, 'pathList') && ~isempty(net.pathList)
    [pIdx, pLoss] = bestPathByHeadLoss(abs(dh), net.pathList);
    results.highestLossPath.pathIndex = pIdx;
    results.highestLossPath.branchIds = net.pathList{pIdx};
    results.highestLossPath.totalLoss = pLoss;
else
    results.highestLossPath = [];
end

fprintf('\nConverged in %d iterations.\n', iterUsed);
if ~isempty(results.Re)
    fprintf('(a) Turbulent branches: %d/%d (Re >= 4000).\n', sum(results.isTurbulent), nPipe);
else
    fprintf('(a) Turbulence check skipped (need D and nu).\n');
end
fprintf('(b) First %d-iteration Qi table stored in results.iterationTableQ.\n', min(nPrint,iterUsed));
fprintf('(c) Converged branch head losses stored in results.dh.\n');
if ~isempty(results.highestLossPath)
    fprintf('(d) Highest-loss path = path #%d, total loss = %.6f.\n', ...
        results.highestLossPath.pathIndex, results.highestLossPath.totalLoss);
else
    fprintf('(d) Highest-loss path skipped (provide net.pathList).\n');
end
end

function K = branchK(net, Q)
if isfield(net, 'K') && ~isempty(net.K)
    K = net.K(:);
    return;
end

required = {'L','D','eps','nu','g'};
for i = 1:numel(required)
    if ~isfield(net, required{i})
        error('Missing net.%s. Provide net.K or full Darcy-Weisbach fields.', required{i});
    end
end

L = net.L(:);
D = net.D(:);
eps = net.eps(:);
nu = net.nu;
g = net.g;

A = pi .* D.^2 ./ 4;
V = abs(Q) ./ A;
Re = max(V .* D ./ nu, 1);

f = zeros(size(Q));
for i = 1:numel(Q)
    if Re(i) < 2300
        f(i) = 64 / Re(i);
    else
        f(i) = 0.25 / (log10(eps(i)/(3.7*D(i)) + 5.74/(Re(i)^0.9))^2);
    end
end

K = (8 .* f .* L) ./ (pi^2 .* g .* D.^5);
end

function out = getFieldOr(s, name, defaultVal)
if isfield(s, name)
    out = s.(name);
else
    out = defaultVal;
end
end

function [bestIdx, bestLoss] = bestPathByHeadLoss(absDh, pathList)
bestIdx = 1;
bestLoss = -inf;
for i = 1:numel(pathList)
    p = pathList{i};
    loss = sum(absDh(p));
    if loss > bestLoss
        bestLoss = loss;
        bestIdx = i;
    end
end
end
