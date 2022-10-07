function obj = calculate_approx_mrpi(obj)
%CALCULATE_MRPI This function generates the minimum RPI set in central form
    obj.n_x = 4;
    obj.n_b = 2;
    obj.delta_blocks = ones(2,2);
    obj.options.solver_sdp = 'SeDuMi';
    q_r_inv = sdpvar(obj.n_x, obj.n_x);
    tau_x = sdpvar();
    tau_w = sdpvar(obj.n_b, 1);
    
    % Matrix representation of S-Procedure variables
    lambda_w = [];
    for i = 1:obj.n_b
        lambda_w = blkdiag(lambda_w, tau_w(i) * eye(obj.delta_blocks(i, 1)));
    end
    
    % Define RPI requirement LMI
    rpi_lmi = blkvar;
    rpi_lmi(1,1) = - q_r_inv;
    rpi_lmi(1,2) = (obj.a - obj.b_u * obj.gcc.k) * q_r_inv;
    rpi_lmi(1,3) = obj.b_w;
    rpi_lmi(2,2) = - tau_x * q_r_inv;
    rpi_lmi(3,3) = - lambda_w;
    
    % Define constraints
    constraints = [
        rpi_lmi <= 0;
        tau_w >= 0;
        tau_x + sum(tau_w) <= 1;
    ];

    % Outer-bounding constraint LMI for each disturbance block
    block_range = [
        0, 0;
        cumsum(obj.delta_blocks, 1);
    ];

    for i = 1:obj.n_b 
        r_y = [block_range(i, 2) + 1, block_range(i + 1, 2)];
        c_y_cl = obj.c_y(r_y(1):r_y(2), :) - obj.d_y_u(r_y(1):r_y(2), :) * obj.gcc.k;
        bound_lmi = blkvar;
        bound_lmi(1, 1) = - eye(size(c_y_cl, 1));
        bound_lmi(1, 2) = c_y_cl * q_r_inv;
        bound_lmi(2, 2) = - q_r_inv;
        constraints = [
            constraints;
            bound_lmi <= 0;
        ];
    end
    
    % Define YALMIP optimization problem
    objective = trace(q_r_inv);
    options = sdpsettings('solver','SeDuMi', 'verbose', 1);
    
    % Solve optimization problem
    opt = optimizer(constraints, objective, options, tau_x, {q_r_inv, tau_w});
    
    tau_list = [];
    cost_list = [];
    q_r_inv_list = [];
    for tau = [0:0.01:1]
        [opt_out, status] = opt(tau);
        
        if status == 0
            tau_list(:, end+1) = [tau, opt_out{2}'];
            cost_list(end+1) = - logdet(opt_out{1});
            if size(q_r_inv_list, 1) == 0
                q_r_inv_list = opt_out{1};
            else
                q_r_inv_list(:, :, end+1) = opt_out{1};
            end
        end
    end
    
    [~, min_idx] = min(cost_list);
    obj.mrpi.q_r_inv = q_r_inv_list(:, :, min_idx);
    obj.mrpi.a_alpha = tau_list(1, min_idx);
    obj.mrpi.a_lambda = tau_list(2:end, min_idx);
end