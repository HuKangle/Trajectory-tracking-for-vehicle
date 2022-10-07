%***************************************************************%
% LinearTireModel_WithKs_FOH
% Deem Vx as constant
%***************************************************************%
function [StateSpaceModel] = Model_Tube_MPC_724(VehiclePara, MPCParameters, Vel, Calpha_f, Calpha_r, CostWeights)
    % generate State-space model
   
M       = VehiclePara.m ;
lf      = VehiclePara.Lr;% 1.12;%1.05%  
lr      = VehiclePara.Lr;% 1.48;%
Iz      = VehiclePara.Iz;
Ts = MPCParameters.Ts;

Theta_1     = Calpha_f + Calpha_r;  
Theta_2     = lf*Calpha_f - lr*Calpha_r;
Theta_3     = lf*lf*Calpha_f + lr*lr*Calpha_r;

Mint = [  M*Vel      0; 
          0          Iz];

Nint = [Theta_1    -M*Vel + Theta_2/Vel;
        Theta_2    Theta_3/Vel];         

F1int = [-Calpha_f; -lf*Calpha_f];

Ac_11     = Mint\Nint; % 2*2
B1cn_11   = Mint\F1int; % 2*1

A_22 =  (Calpha_f + Calpha_r)/ Vel/ M;
A_23 = -(Calpha_f +Calpha_r)/ M;
A_24 =  (lf * Calpha_f - lr * Calpha_r)/ Vel/ M;
A_42 = (lf * Calpha_f - lr * Calpha_r)/ Vel/ Iz;
A_43 = -(lf * Calpha_f - lr * Calpha_r)/ Iz;
A_44 = (lf^2 * Calpha_f + lr^2 * Calpha_r)/ Vel/ Iz;

ac_12   = zeros(2,4);
ac_21   = zeros(4,2);     
  
ac_22   =  [0 1 0 0;
            0 A_22 A_23 A_24;
            0 0 0 1;
            0 A_42 A_43 A_44];
acn      = [Ac_11 ac_12; ac_21 ac_22];

b1cn_13 = [0; -Calpha_f/M; 0; -lf*Calpha_f/Iz];
b2cn_13 = [0; A_24*Vel - Vel^2; 0; A_44 * Vel];

b1cn    = [B1cn_11; b1cn_13];
b2cn    = [zeros(2,1); b2cn_13];

N_x = size(acn,2);
N_u = size(b1cn,2);
N_w = size(b2cn,2);
Mc = [acn, b1cn, b2cn;
      zeros(N_u + N_w, N_x + N_u + N_w)];
Md = expm(Ts * Mc);
acd1 = Md(1:N_x, 1:N_x);
bcd1 = Md(1:N_x, N_x + N_u); 
bcd2 = Md(1:N_x, N_x + N_u + 1 : N_x + N_u + N_w);

%% lqr for error-based model
% states as [ed_dot, ed, epsi_dot, epsi]
% input as delta_f
% states obtained by the vehicle position and the planned position
% Calpha_f = VehiclePara.CafHat;
% Calpha_r = VehiclePara.CarHat;
A_22 = (Calpha_f + Calpha_r)/ Vel/ M;
A_23 = -(Calpha_f +Calpha_r)/ M;
A_24 = (lf * Calpha_f - lr * Calpha_r)/ Vel/ M;
A_42 = (lf * Calpha_f - lr * Calpha_r)/ Vel/ Iz;
A_43 = -(lf * Calpha_f - lr * Calpha_r)/ Iz;
A_44 = (lf^2 * Calpha_f + lr^2 * Calpha_r)/ Vel/ Iz;
A_e = [0 1 0 0;
        0 A_22 A_23 A_24;
        0 0 0 1;
        0 A_42 A_43 A_44];
B_e = [0; -Calpha_f/M; 0; -lf*Calpha_f/Iz];

%% State-Space Model for error-based model
a = [0 0 1 0; 
     0 0 0 1;
     0, A_23, A_22, A_24;
     0, A_43, A_42, A_44];
b_u = [0; 0;  -Calpha_f/M; -lf*Calpha_f/Iz];
b_r = [0; 0; (lf * Calpha_f - lr * Calpha_r)/ M - Vel^2; (lf^2 * Calpha_f + lr^2 * Calpha_r)/ Iz];

dC = [Calpha_f Calpha_r]*0.3;
dCf = dC(1);
dCr = dC(2);
In   = diag([1/M, 1/Iz]);
LArm = [1,  1;
        lf, -lr];
Larm = [-1,1;
        1, 1;
        lf, -lr ];
DMag = diag([dCf, dCr]);
S    = [1; 0];

H  = In * LArm;
Ea = - DMag * Larm' / Vel;
Eb = DMag * S;
Ec = DMag * [lf; -lr];

H = [zeros(2); H];
Ea = [zeros(2,1) Ea];
Ebu = Eb;
Ebr = Ec;

Nx = size(a, 2);
Nu = size(b_u, 2);
Nr = size(b_r,2);
Nw = size(H, 2);
N_x = Nx;
N_u = Nu;
N_r = Nr;
N_w = Nw;
Mc = [a, b_u, b_r, H;
      zeros(N_u + N_r+ N_w, N_x + N_u + N_r + N_w)];
Md = expm(Ts * Mc);
F = Md(1:N_x, 1:N_x);
G = Md(1:N_x, N_x + N_u); 
b_r = Md(1:N_x, N_x + N_u + 1 : N_x + N_u + N_r);
b_w = Md(1:N_x, N_x + N_u + N_r + 1 : N_x + N_u + N_r + N_w);
 Ef = Ea * Ts;
 Eg = Eb * Ts;
 Ebr = Ec * Ts;
 Q = diag([CostWeights.Wey CostWeights.Wephi CostWeights.Weyd CostWeights.Wepsid]);
 R  = CostWeights.Sl;
 b_w = b_w*1000;
 Ef  = Ef/1000;
 Eg  = Eg/1000;
 Ebr = Ebr /1000;
 kPosDefTest = 1e-10;   % Minimum eigenvalue to consider matrix Positive-Definit
  n =zeros(Nx, Nu);
    w = [Q n;
        n' R];
    % Decompose weight matrix
    [~, s, vt] = svd(w);
    mask = (diag(s) >= kPosDefTest);
    w_half = sqrt(s(mask, mask)) * vt(:, mask)';
	c_z = w_half(:,1:Nx);
    d_z_u = w_half(:,Nx+1:end);   
    % Set instance variables
%% ZOH/FOH
    StateSpaceModel.AC      = Ac_11;
    StateSpaceModel.BC1     = B1cn_11;
    StateSpaceModel.acn     = acn;
    StateSpaceModel.b1cn    = b1cn; 
    StateSpaceModel.b2cn    = b2cn;  
    
    StateSpaceModel.acd1     = acd1;
    StateSpaceModel.bcd1     = bcd1; 
    StateSpaceModel.bcd2     = bcd2; 

    % new state space model
    StateSpaceModel.A_e     = A_e;
    StateSpaceModel.B_e     = B_e;
    
    % Error-based Model
    StateSpaceModel.a       = F;
    StateSpaceModel.b_u     = G;
    StateSpaceModel.b_r     = b_r;
    StateSpaceModel.b_w     = b_w;
    StateSpaceModel.c_y     = Ef;
    StateSpaceModel.d_y_u   = Eg;
    StateSpaceModel.d_w_u   = Ebr;

    StateSpaceModel.c_z     = c_z;
    StateSpaceModel.d_z_u   = d_z_u;
    StateSpaceModel.n_z     = size(c_z, 1);
    StateSpaceModel.n_w     = Nw;
    StateSpaceModel.n_x     = Nx;
    StateSpaceModel.n_u     = Nu;
    StateSpaceModel.n_b     = Nu;

    StateSpaceModel = gcc(StateSpaceModel);
end % end of func_SpatialDynamicalModel

 