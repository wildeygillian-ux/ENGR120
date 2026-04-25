% 22-branch assignment runner using the pipe properties provided in the prompt.
% Fill in Q0, C (loop matrix), and pathList according to your network diagram.

clear; clc;

% ----------------------------
% Given branch data (from prompt)
% ----------------------------
net.L = [ ...
    250 100 200 200 250 500 150 100 100 100 150 100 ...
    200 150 100 150  50 200 250 200 250 200]';          % m

net.D = [ ...
    15 10 10 10 10 15 10 10 5 5 5 5 10 5 5 5 5 10 10 10 10 10]' ./ 100;  % m

net.eps = [ ...
    0.2 0.15 0.15 0.15 0.2 0.2 0.2 0.2 0.1 0.1 0.1 0.1 ...
    0.15 0.1 0.1 0.1 0.1 0.15 0.2 0.2 0.1 0.2]' ./ 1000;                  % m

% Fluid properties (water)
net.nu = 8.93e-7;      % m^2/s
net.g = 9.81;          % m/s^2

% -----------------------------------------------
% REQUIRED from your network diagram / assignment
% -----------------------------------------------
% 1) Initial flow guess in each branch (22x1)
net.Q0 = [ ...
    0.05 0.04 0.03 0.02 0.02 0.03 0.01 0.01 0.005 0.005 0.005 ...
    0.005 0.01 0.005 0.005 0.005 0.004 0.01 0.02 0.01 0.01 0.01]';

% 2) Loop-incidence matrix C (nLoops x 22), +1/-1/0 by loop direction
%    Replace the placeholder below with the actual loop definitions.
net.C = zeros(1,22);
net.C(1,[1 2 3]) = [1 1 -1];   % placeholder only

% 3) Candidate source-to-demand paths using branch IDs (for part d)
%    Replace with actual branch sequences from your diagram.
net.pathList = {
    [1 2 3], ...
    [4 5 6]
};

net.maxIter = 200;
net.tol = 1e-10;
net.nPrint = 4;

% Solve
results = solve_pipe_network(net);

% --------------------
% (a) Turbulence check
% --------------------
if isempty(results.Re)
    error('Reynolds number output missing. Ensure net.D and net.nu are set.');
end
TRe = table((1:22)', results.Re, results.isTurbulent, ...
    'VariableNames', {'Pipe','Re','Turbulent'});
disp('Part (a): Reynolds/turbulence check by branch');
disp(TRe);

% -------------------------------
% (b) First 4 iterations, Qi table
% -------------------------------
rowNames = compose('Iter%d', 1:size(results.iterationTableQ,1));
varNames = compose('Q%d', 1:22);
TQ = array2table(results.iterationTableQ, 'RowNames', rowNames, 'VariableNames', varNames);
disp('Part (b): First 4 iterations of Qi (Table-3 style)');
disp(TQ);

% -----------------------------------
% (c) Converged branch head losses dh
% -----------------------------------
Tdh = table((1:22)', results.Q, results.dh, ...
    'VariableNames', {'Pipe','Q','Delta_h'});
disp('Part (c): Converged flow and head loss in each branch');
disp(Tdh);

% ------------------------------
% (d) Highest head-loss path
% ------------------------------
if isempty(results.highestLossPath)
    warning('No pathList provided for part (d).');
else
    fprintf('Part (d): Highest head-loss path index = %d\n', results.highestLossPath.pathIndex);
    fprintf('Branch IDs: %s\n', mat2str(results.highestLossPath.branchIds));
    fprintf('Total path head loss: %.6f m\n', results.highestLossPath.totalLoss);
end
