function [RMPC_controller] = Solver_HTMPC(MPCParameters, Constraints, CostWeights, obj)
%=====================================%
ts     = MPCParameters.Ts;
Rdun	 = CostWeights.Rdun;  

% obj = StateSpaceModel;
dumax    = Constraints.dumax;
umax     = Constraints.umax;  
DPhimax = Constraints.DPhimax;
Dymax   = Constraints.Dymax;
h_x      = [1 0 0 0 ;
            -1 0 0 0; 
            0 1 0 0 ;
            0 -1 0 0; 
            0 0 1 0 ;
            0 0 -1 0;
            0 0 0  1;
            0 0 0 -1; 
            zeros(2,4) ];
h_u      = [zeros(8,1); 1; -1];
g = - [Dymax ; Dymax; DPhimax; DPhimax; 3; 3; 3; 3; umax; umax];
a_alpha  = obj.mrpi.a_alpha;
a_lambda = obj.mrpi.a_lambda;
q_r_inv  = obj.mrpi.q_r_inv;

obj.n_c     = size(h_x,1);
obj.n_t_c   = obj.n_c;
obj.h_t_x   = h_x;
obj.g_t     = g;
kSlackWeight = 500;
obj.kSlackWeight  = kSlackWeight;
obj.kAlphaWeight =100; 
% Set horizon
obj.n_t = MPCParameters.Np;
% obj.n_t = 5;

% Enumarate vertexes of A and B
a = zeros(size(obj.a));
b_u = zeros(size(obj.b_u));
b_r = zeros(size(obj.b_u));
disturbance_edges = delta_vector(obj.n_w);
n_systems = size(disturbance_edges, 1);
for i = 1 : n_systems
    delta = diag(disturbance_edges(i, :));
    a(:, :, i) = obj.a + obj.b_w * delta * obj.c_y;
    b_u(:, :, i) = obj.b_u + obj.b_w * delta * obj.d_y_u;
    b_r(:, :, i) = obj.b_r + obj.b_w * delta * obj.d_w_u;
end

% Define optimization variables
z = sdpvar(obj.n_x, obj.n_t + 1, 'full');
v = sdpvar(obj.n_u, obj.n_t, 'full');
w_c = sdpvar(obj.n_u, obj.n_t, 'full');
alpha = sdpvar(1, obj.n_t + 1, 'full');
du = sdpvar(obj.n_u, obj.n_t, 'full');
z0 = sdpvar(obj.n_x, 1, 'full');
pastu = sdpvar(obj.n_u, 1, 'full');
% Define problem objective
p = obj.gcc.p;
w_half = [obj.c_z, obj.d_z_u];
w = w_half' * w_half;
Objective = 0;
for i = 1:obj.n_t
    Objective = Objective + [z(:,i); v(:, i)]' * w * [z(:, i); v(:,i)] + obj.kAlphaWeight * (0.5-alpha(i))^2;
end
Objective = Objective + z(:, obj.n_t + 1)' * p * z(:, obj.n_t + 1) + obj.kAlphaWeight *(0.5-alpha(obj.n_t+1))^2;

% Define system dynamics constraints
constraint = [];
e_r_half = q_r_inv ^ -0.5;
% a_alpha = 0.9;
for i = 1 : n_systems
    for j = 1 : obj.n_t
        a_tilda_i = a(:, :, i) - b_u(:, :, i) * obj.gcc.k;
        a_alpha = norm(e_r_half * a_tilda_i / e_r_half, 2);
        dynamics_error = a(:, :, i) * z(:,j) + b_u(:, :, i) * v(:, j) + b_r(:,:,i) * w_c(:,j) - z(:,j + 1);
        constraint = [
            constraint;
            a_alpha * alpha(j) + norm(e_r_half * dynamics_error, 2) <= alpha(j + 1);
        ];
    end
end
constraint = [
    constraint;
    alpha >= 0;
];
h_tilda = h_x - h_u * obj.gcc.k;
h_alpha = sqrt(sum((h_tilda * (q_r_inv ^ 0.5)) .^ 2, 2));
h_value = h_x * z(:,1:end-1) + h_u * (v ) + g * ones(1,  obj.n_t);
slack = sdpvar(obj.n_c, obj.n_t);
    
% Objective
Objective = Objective + obj.kSlackWeight * sum(sum(slack));

% Constraint slacking
constraint = [
    constraint;
    h_value + h_alpha * alpha(1:end-1) <= slack;
    slack >= 0;];
% terminal constraint
h_t_alpha = sqrt(sum((obj.h_t_x * (q_r_inv ^ 0.5)) .^ 2, 2));
h_t_value = obj.h_t_x * z(:,end) + obj.g_t;
slack = sdpvar(obj.n_t_c, 1);
% Objective
Objective = Objective + obj.kSlackWeight * sum(slack);
constraint = [
    constraint;
    h_t_value + h_t_alpha * alpha(end) <= slack;
    slack >= 0;
];

for i = 1 : obj.n_t
     if i == 1
      constraint = [ constraint;
                    ts * du(:,i) == v(:,i) - pastu];
    else
      constraint = [constraint;
                    ts * du(:,i) == v(:,i) - v(:,i-1)];
     end
     constraint = [constraint ; -dumax <= du(:,i) <= dumax;];
end

Objective = Objective + du * Rdun * du';

parameters_in = {z(:,1), w_c, pastu,alpha(1)};
outputs = {z, v, alpha, Objective};
ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
% problem = optimize(constraint, Objective, ops);
RMPC_controller = optimizer(constraint, Objective, ops, parameters_in, outputs);
end
% x_opt = double(z);
% u_opt = double(v);
% objective = double(Objective);
% controller = optimizer(constraint, Objective, ops, z(:,1), outputs);
% end
function deltas = delta_vector(k)
    values = [-1, 1];
    
    n = numel(values);
    combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
    combs = reshape(values(combs),[],k);
    
    deltas = [];
    for i = 1 : size(combs, 1)
        deltas = [deltas; unique(perms(combs(i, :)), 'rows')];
    end
end