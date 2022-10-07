function [RMPC_controller] = Solver_TMPC(VehiclePara, MPCParameters, Constraints, CostWeights, StateSpaceModel, Envelope, Vel)
    % YalmipInterface

    np = MPCParameters.Np;
%     np= 5;
%     nc = 5;
    nx = 4;
    nu = 1;
    As       =  StateSpaceModel.Ad1;
    Bs1      =  StateSpaceModel.Bd1;
    Bs2      =  StateSpaceModel.Bd2;
    dumax    = Constraints.dumax;
    umax     = Constraints.umax; 
    %=====================================%
    ts       = MPCParameters.Ts;
    Sl       = CostWeights.Sl;
    Rdun	 = CostWeights.Rdun;  
    Wshl	 = 1000; 
    Wend     = 1000;
    Wenv     = 1000;

    Robust_constraint = Constraints.Robust_constraint.TMPC;
    X_bar1      = Robust_constraint.X_bar1;
    X_bar2      = Robust_constraint.X_bar2;
    X_f         = Robust_constraint.X_f;
    Q           = CostWeights.Q_err;
    P           = Robust_constraint.terminal_cost;
    N1       = size(X_bar1.b,1);
    N2       = size(X_bar2.b,1);
    %% Building the object controller
    x0 = sdpvar(nx, 1);
    x_nom    = sdpvar(nx, np+1);
    u        = sdpvar(nu,  np);
    du       = sdpvar(nu,  np);
    pastu    = sdpvar(nu, 1);
    Uaug     = sdpvar(2, np);
    S_env    = sdpvar(N1, np);
    Sh       = sdpvar(N2, np);
    S_term   = sdpvar(1, 1);

    objective = 0;

    for i = 1: np
        if i == 1
          constraints = [ x_nom(:,1) == x0;
                         ts * du(:,i) == u(:,i) - pastu];
        else
            constraints = [constraints;ts * du(:,i) == u(:,i) - u(:,i-1)];
        end
        
        constraints = [ constraints;
        x_nom(:,i+1) == As * x_nom(:,i) + Bs1* u(:,i) - Bs2 * Uaug(2,i);];
    
        constraints = [ constraints; 
            -umax <= u(:,i) <= umax; 
            -dumax <= du(:,i) <= dumax; 
             X_bar1.A * x_nom(1:2,i) <= X_bar1.b +  S_env(:,i);
             X_bar2.A * x_nom(3:4,i) <= X_bar2.b + Sh(:,i);
             0 <= S_env(:,i);
             0 <= Sh(:,i);
             ];
    end
    % terminal constraints
    constraints = [constraints; 
                   X_f.A * x_nom(3:4,np+1) <= X_f.b + S_term;
                   0 <= S_term;
                   ];
    % objective function
    Q = diag([0.1,0.1,75,50]);
    for i = 1 : np
    objective = objective + du(:,i)'*Rdun*du(:,i) + u(:,i)'*Sl*u(:,i);
    objective = objective + Sh(:,i)'* Wshl *Sh(:,i) + S_env(:,i)'*Wenv*S_env(:,i);
    objective = objective + x_nom(:,i)' * Q * x_nom(:,i) ;
    end
    % terminal cost
  objective = objective + x_nom(:,np+1)' * P * x_nom(:,np+1);
% Solver options
    ops = sdpsettings('solver','sedumi','verbose',0);
    parameters_in = {x0, Uaug, pastu};
    solutions_out = {x_nom,u,du, objective};
    RMPC_controller = optimizer(constraints, objective, ops, parameters_in, solutions_out);
