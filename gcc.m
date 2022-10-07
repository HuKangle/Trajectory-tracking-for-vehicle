
function obj = gcc(obj)
%CALCULATE_GCC This function calculates the linear Guaranteed Cost Controller for the nominal region

    obj.n_y = size(obj.c_y, 1);
    obj.delta_blocks = [obj.n_w , obj.n_y];
    
    % Define LMI variables
    p_inv = sdpvar(obj.n_x, obj.n_x);            % P^(-1)
    k_p_inv = sdpvar(obj.n_u, obj.n_x, 'full');  % K P^(-1)
    s = sdpvar(obj.n_x, obj.n_x);                % S >= P
    e = sdpvar(obj.n_b, 1);
    % Matrix representation of S-Procedure variables
    upsilon_w = [];
    upsilon_y = [];
    for i = 1:obj.n_b
        upsilon_w = blkdiag(upsilon_w, e(i) * eye(obj.delta_blocks(i, 1)));
        upsilon_y = blkdiag(upsilon_y, e(i) * eye(obj.delta_blocks(i, 2)));
    end
    
    % Define GCC robustness requirement LMI
    gcc_lmi = blkvar;
    gcc_lmi(1,1) = -upsilon_y;
    gcc_lmi(1,4) = obj.c_y * p_inv - obj.d_y_u * k_p_inv;
    gcc_lmi(2,2) = -eye(obj.n_z);
    gcc_lmi(2,4) = obj.c_z * p_inv - obj.d_z_u * k_p_inv;
    gcc_lmi(3,3) = - p_inv + obj.b_w * upsilon_w * obj.b_w';
    gcc_lmi(3,4) = obj.a * p_inv - obj.b_u * k_p_inv;
    gcc_lmi(4,4) = - p_inv;
    
    % Define GCC cost LMI
    cost_lmi = blkvar;
    cost_lmi(1,1) = - s;
    cost_lmi(1,2) = eye(obj.n_x);
    cost_lmi(2,2) = - p_inv;
    
    % Define YALMIP optimization problem
    constraints = [gcc_lmi <= 0; 
                   cost_lmi <= 0];
    objective = trace(s);
    options = sdpsettings('solver', 'sedumi', 'verbose', 0);
    
    % Solve optimization problem
    solve_out = optimize(constraints, objective, options);
    
    if solve_out.problem ~= 0
        error('SDP solver did not converge, please check if your problem is correct');
    end
    
    % S-Procedure variable
    upsilon_w = value(upsilon_w);
    upsilon_y = value(upsilon_y);
    
    % GCC cost matrix
    obj.gcc.p = inv(value(p_inv));
    % GCC gain matrix
    obj.gcc.k = value(k_p_inv) * obj.gcc.p;
    % Suboptimal gains, needed to get r_bar
    obj.gcc.x = inv(value(p_inv) - obj.b_w * upsilon_w * obj.b_w');
    % For a controller u = - K x + v, the new cost matrix is v' * r_bar * v
    obj.gcc.r_bar = (obj.d_z_u' * obj.d_z_u) + ...
                    obj.d_y_u' * upsilon_y * obj.d_y_u + ...
                    obj.b_u' * obj.gcc.x * obj.b_u;
    % Set calculation flag
    obj.is_gcc_set = true;
end

