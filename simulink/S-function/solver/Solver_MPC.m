function [RMPC_controller] = Solver_MPC(VehiclePara, MPCParameters, Constraints, CostWeights, StateSpaceModel, Envelope, Vel)
    % YalmipInterface
    %===============================%  MPC parameters and model 
    lf = VehiclePara.Lf;
    lr = VehiclePara.Lr;
    M  = VehiclePara.m;
    g  = VehiclePara.g;

    CarHat = -92000;
%     nx = MPCParameters.Nx;
%     nu = MPCParameters.Nu;
    np = MPCParameters.Np;
%     ns = MPCParameters.Ns;

    nx = 8;
    nu = 1;
    
    As       =  StateSpaceModel.acd1;
    Bs1      =  StateSpaceModel.bcd1;
    Bs2      =  StateSpaceModel.bcd2;

    %=====================================%
    ts     = MPCParameters.Ts;
    tl     = MPCParameters.Tsl;
    Q        = CostWeights.Q_err;
    Sl       = CostWeights.Sl;
    Rdun	 = CostWeights.Rdun;  
%     Rdul	 = CostWeights.Rdul;  
%     Wshl	 = CostWeights.Wshl; 
%     Wshr     = CostWeights.Wshr;
%     Wenv     = diag([Wshr,Wshr]);
%     Yshl     = 1000;
%     Sshmax   = Constraints.Sshmax;
    %=======================================%
    dumax    = Constraints.dumax;
    umax     = Constraints.umax;  
    Alphar_lim = Constraints.arlim;

%     Np          = MPCParameters.Np;

    Henv  = [1 ; -1]; 
    Hsh  = [1/Vel   -lr/Vel    0   0    ;
            0         1        0   0    ];   
%     Psh  = [0  0;  g/Vel   0];
    
    r_ssmax = -CarHat*Alphar_lim*(1+lr/lf)/(M*Vel);                  

    Gsh  = [Alphar_lim    r_ssmax]';   
    Genv  = Envelope.Genv;
%     Hyzmp    = Envelope.H_yzmp_4;
%     Pyzmp1   = Envelope.P_yzmp1_4;
%     Pyzmp2   = Envelope.P_yzmp2_4;

%     YZPM_lim = 0.8;
%     LTR_lim  = 0.9;
%     K_phi    = 145330;  % roll stiffness coefficient
%     D_phi    = 4500;    % roll damping coefficient
%     H_LTR    = 2 / (VehiclePara.m*VehiclePara.g * VehiclePara.Tr) * [D_phi, K_phi];
    % H3       = H_LTR(3);
    % H4       = H_LTR(4);

    %======================================%
    %% Building the object controller
    x0 = sdpvar(nx, 1);
    x_nom  = sdpvar(nx, np+1);
    u = sdpvar(nu,  np);
    du = sdpvar(nu,  np);
    pastu = sdpvar(nu, 1);
    Uaug = sdpvar(2, np);
%     Sh   = sdpvar(2, np);
%     Yh   = sdpvar(1, np);
%     S_env = sdpvar(2, np);
%     F_yf  = sdpvar(1,np);
    objective = 0;
    %  STATE: [beta, yawrate, DPHI, PHI, ey, ey_dot, epsi, epsi_dot]
    %  First phase setting
    for i = 1: np
        if i == 1
          constraints = [ x_nom(:,1) == x0;
                         ts * du(:,i) == u(:,i) - pastu];
        else
            constraints = [constraints;ts * du(:,i) == u(:,i) - u(:,i-1)];
        end
        
        constraints = [ constraints;
        x_nom(:,i+1) == As * x_nom(:,i) + Bs1* u(:,i) + Bs2 * Uaug(:,i);];
    
        constraints = [ constraints; 
            -umax <= u(:,i) <= umax; 
            -dumax <= du(:,i) <= dumax; 
%              Henv * x_nom(5,i) <= Genv{i,1} + S_env(:,i);
              Henv * x_nom(5,i) <= Genv{i,1} ;
%             -Gsh - Sh(:,i)  <= Hsh * x_nom(1:4,i) <= Gsh + Sh(:,i)] ; 
             -Gsh  <= Hsh * x_nom(1:4,i) <= Gsh;] ; 
    end
    % terminal constraints
    % constraints = [ constraints; -LTR_lim <= H_LTR * x_nom(3:4,np) <= LTR_lim];    
    %% objective function
    Q = diag([1000,50,50,100]);
    for i = 1 : np
    objective = objective + du(:,i)'*Rdun*du(:,i) + u(:,i)'*Sl*u(:,i);
%     objective = objective + Sh(:,i)'* Wshl *Sh(:,i) + S_env(:,i)'*Wenv*S_env(:,i);
    objective = objective + x_nom(5:8,i)' * Q * x_nom(5:8,i) ;
    end
% Solver options
    ops = sdpsettings('solver','gurobi','verbose',0);
    parameters_in = {x0, Uaug, pastu};
    solutions_out = {x_nom,u,du, objective};
    RMPC_controller = optimizer(constraints, objective, ops, parameters_in, solutions_out);
