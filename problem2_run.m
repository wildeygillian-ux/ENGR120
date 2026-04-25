% Problem 2 - 22-branch community pipe network (Figure 2)
% This script is fully parameterized with the provided L, d, epsilon data
% and the given inlet/outlet flow boundary conditions.

clear; clc;

data = problem2_data();
results = solve_problem2_network(data);

% (a) Turbulence verification
Ta = table((1:22)', results.Re, results.isTurbulent, ...
    'VariableNames', {'Pipe','Re','Turbulent'});
disp('Part (a): Turbulence check in each branch (Re >= 4000)');
disp(Ta);

% (b) First 4 iterations of Qi (Table-3 style)
Tb = array2table(results.iterationTableQ, ...
    'RowNames', compose('Iter%d', 1:size(results.iterationTableQ,1)), ...
    'VariableNames', compose('Q%d', 1:22));
disp('Part (b): First 4 iterations of branch flow rates Qi [m^3/s]');
disp(Tb);

% (c) Converged head loss per branch
Tc = table((1:22)', results.Q, results.dh, ...
    'VariableNames', {'Pipe','Q_m3s','Delta_h_m'});
disp('Part (c): Converged branch flow and head loss');
disp(Tc);

% (d) Highest head-loss source-to-outlet path
fprintf('Part (d): Highest head-loss outlet path = %.6f m\n', results.highestLossPath.totalLoss);
fprintf('Outlet node = %d\n', results.highestLossPath.outletNode);
fprintf('Path nodes = %s\n', mat2str(results.highestLossPath.nodes));
fprintf('Path branches = %s\n', mat2str(results.highestLossPath.branches));

function data = problem2_data()
% Pipe properties from Problem2 Pipes.xls (as provided)
data.L = [250 100 200 200 250 500 150 100 100 100 150 100 ...
          200 150 100 150 50 200 250 200 250 200]';
data.D = [15 10 10 10 10 15 10 10 5 5 5 5 10 5 5 5 5 10 10 10 10 10]' / 100; % m
data.eps = [0.2 0.15 0.15 0.15 0.2 0.2 0.2 0.2 0.1 0.1 0.1 0.1 ...
            0.15 0.1 0.1 0.1 0.1 0.15 0.2 0.2 0.1 0.2]' / 1000; % m

% Fluid properties
data.nu = 8.93e-7; % m^2/s
data.g = 9.81;     % m/s^2

% Network node indexing for Fig. 2 interpretation (16 nodes)
% Branch i connects node from(i) -> to(i)
data.from = [1 2 7 8 6 1 2 3 4 10 7 3 10 8 11 12 13 14 9 4 5 13]';
data.to   = [2 7 8 9 9 6 3 4 12 12 10 10 11 11 14 13 14 15 15 5 16 16]';

% External nodal flows (+ demand/outlet, - supply/inlet)
data.nNodes = 16;
data.demand = zeros(data.nNodes,1);
data.demand(1)  = -0.10; % inlet Qin
data.demand(11) =  0.02; % outlet 1
data.demand(16) =  0.03; % outlet 2
data.demand(15) =  0.05; % outlet 3

% Solver controls
data.refNode = 1;
data.refHead = 0;
data.head0 = zeros(data.nNodes,1);
data.maxIter = 100;
data.tol = 1e-10;
data.qTol = 1e-12;
data.qMaxIter = 100;

% Source and outlets for path analysis
data.sourceNode = 1;
data.outletNodes = [11 15 16];
end
