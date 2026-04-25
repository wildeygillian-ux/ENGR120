function results = solve_problem2_network(data)
%SOLVE_PROBLEM2_NETWORK Nodal-head solver for Problem 2 (22-branch network).

from = data.from(:);
to = data.to(:);
L = data.L(:);
D = data.D(:);
eps = data.eps(:);
nu = data.nu;
g = data.g;

nBranch = numel(from);
nNodes = data.nNodes;
if nBranch ~= 22
    warning('Expected 22 branches, got %d.', nBranch);
end

head = data.head0(:);
head(data.refNode) = data.refHead;
unknown = setdiff((1:nNodes)', data.refNode);

maxIter = data.maxIter;
tol = data.tol;
qTol = data.qTol;
qMaxIter = data.qMaxIter;

Qhist = zeros(min(4,maxIter), nBranch);

for k = 1:maxIter
    [Q, ~, Re, dh] = branchFlowsFromHeads(head, from, to, L, D, eps, nu, g, qTol, qMaxIter);

    if k <= 4
        Qhist(k,:) = Q.';
    end

    R = continuityResidual(Q, from, to, data.demand, nNodes);
    Rred = R(unknown);

    if norm(Rred, inf) < tol
        break;
    end

    J = zeros(numel(unknown));
    hStep = 1e-6;
    for j = 1:numel(unknown)
        hp = head;
        hp(unknown(j)) = hp(unknown(j)) + hStep;
        [Qp,~,~,~] = branchFlowsFromHeads(hp, from, to, L, D, eps, nu, g, qTol, qMaxIter);
        Rp = continuityResidual(Qp, from, to, data.demand, nNodes);
        J(:,j) = (Rp(unknown)-Rred)/hStep;
    end

    dH = -J\Rred;
    head(unknown) = head(unknown) + dH;
    head(data.refNode) = data.refHead;

    if norm(dH,inf) < tol
        break;
    end
end

[Q, f, Re, dh] = branchFlowsFromHeads(head, from, to, L, D, eps, nu, g, qTol, qMaxIter);

results.iterations = k;
results.iterationTableQ = Qhist(1:min(4,k),:);
results.Q = Q;
results.f = f;
results.Re = Re;
results.isTurbulent = Re >= 4000;
results.dh = dh;
results.head = head;

results.highestLossPath = highestLossDirectedPath(data.sourceNode, data.outletNodes, from, to, Q, dh);

fprintf('\nConverged in %d iterations.\n', k);
fprintf('(a) Turbulent branches: %d/%d\n', sum(results.isTurbulent), nBranch);
fprintf('(b) First %d iteration rows stored in results.iterationTableQ\n', size(results.iterationTableQ,1));
fprintf('(c) Converged branch losses stored in results.dh\n');
fprintf('(d) Highest-loss source-to-outlet path: %.6f m\n', results.highestLossPath.totalLoss);
end

function [Q, f, Re, dh] = branchFlowsFromHeads(head, from, to, L, D, eps, nu, g, qTol, qMaxIter)
n = numel(from);
Q = zeros(n,1); f = zeros(n,1); Re = zeros(n,1);
dh = head(from)-head(to);
A = pi.*D.^2/4;
for i=1:n
    s = sign(dh(i));
    if s == 0
        continue;
    end
    dhAbs = abs(dh(i));
    R0 = 8*0.02*L(i)/(g*pi^2*D(i)^5);
    q = s*sqrt(max(dhAbs/R0,0));
    for it=1:qMaxIter
        v = abs(q)/A(i);
        Rei = max(v*D(i)/nu,1);
        fi = frictionFactor(Rei, eps(i), D(i));
        Ri = 8*fi*L(i)/(g*pi^2*D(i)^5);
        qNew = s*sqrt(max(dhAbs/Ri,0));
        if abs(qNew-q) < qTol
            q = qNew;
            break;
        end
        q = qNew;
    end
    Q(i)=q;
    v = abs(q)/A(i);
    Re(i)=v*D(i)/nu;
    f(i)=frictionFactor(max(Re(i),1), eps(i), D(i));
end
end

function f = frictionFactor(Re, eps, D)
if Re < 2300
    f = 64/Re;
else
    f = 0.25/(log10(eps/(3.7*D) + 5.74/(Re^0.9))^2);
end
end

function R = continuityResidual(Q, from, to, demand, nNodes)
R = -demand(:);
for i=1:numel(Q)
    R(from(i)) = R(from(i)) - Q(i);
    R(to(i)) = R(to(i)) + Q(i);
end
R = R(1:nNodes);
end

function out = highestLossDirectedPath(sourceNode, outletNodes, from, to, Q, dh)
nNodes = max(max(from), max(to));
adj = cell(nNodes,1);
edge = cell(nNodes,1);
wt = cell(nNodes,1);
for i=1:numel(Q)
    if Q(i) >= 0
        u = from(i); v = to(i);
    else
        u = to(i); v = from(i);
    end
    adj{u}(end+1) = v; %#ok<AGROW>
    edge{u}(end+1) = i; %#ok<AGROW>
    wt{u}(end+1) = abs(dh(i)); %#ok<AGROW>
end

bestLoss = -inf;
bestNodes = [];
bestBranches = [];
bestOutlet = outletNodes(1);

for t = outletNodes(:)'
    [nodes, branches, loss] = dfsMaxLoss(sourceNode, t, adj, edge, wt);
    if ~isempty(nodes) && loss > bestLoss
        bestLoss = loss;
        bestNodes = nodes;
        bestBranches = branches;
        bestOutlet = t;
    end
end

out.totalLoss = bestLoss;
out.nodes = bestNodes;
out.branches = bestBranches;
out.outletNode = bestOutlet;
end

function [bestNodes, bestBranches, bestLoss] = dfsMaxLoss(src, dst, adj, edge, wt)
visited = false(numel(adj),1);
bestLoss = -inf; bestNodes = []; bestBranches = [];

    function dfs(u, nodes, branches, loss)
        visited(u) = true;
        nodes = [nodes, u]; %#ok<AGROW>
        if u == dst
            if loss > bestLoss
                bestLoss = loss;
                bestNodes = nodes;
                bestBranches = branches;
            end
            visited(u) = false;
            return;
        end
        for k = 1:numel(adj{u})
            v = adj{u}(k);
            if ~visited(v)
                dfs(v, nodes, [branches, edge{u}(k)], loss + wt{u}(k));
            end
        end
        visited(u) = false;
    end

dfs(src, [], [], 0);
end
