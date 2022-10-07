function [h_W] = support_function(a, W)
% by Lukas Brunke, lukas.brunke@mail.utoronto.ca
%Function h_W(a) from eq. (9) from Rakovic et al. (2005).
w = sdpvar(W.Dim, 1);

% Set cost to negative in order to maximize
cost = - a' * w;
constraints = W.A * w <= W.b;

solver = 'linprog';
% solver = 'gurobi';
options = sdpsettings('verbose',0,'solver',solver);
options.OptimalityTolerance = 1e-15;
options.StepTolerance = 1e-15;

% Solve optimization problem
problem = optimize(constraints, cost, options);

% Reverse the sign from the result in order to get the real value 
h_W = - double(cost);
end

