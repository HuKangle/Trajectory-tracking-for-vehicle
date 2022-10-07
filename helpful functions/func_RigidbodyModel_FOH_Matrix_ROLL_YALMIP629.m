%***************************************************************%
% LinearTireModel_WithKs_FOH
% Deem Vx as constant
%***************************************************************%
function [StateSpaceModel] = func_RigidbodyModel_FOH_Matrix_ROLL_YALMIP629(VehiclePara, MPCParameters, Vel, Calpha_f, Calpha_r)
    % generate State-space model
   
M       = VehiclePara.m ;
g       = VehiclePara.g;%m/s^2
L       = VehiclePara.L;  %a = 1.11;  
lf      = VehiclePara.Lr;% 1.12;%1.05%  
lr      = VehiclePara.Lr;% 1.48;%
Tr      = VehiclePara.Tr ; %m
hCG     = VehiclePara.hCG;%m
Ix      = VehiclePara.Ix;
Iz      = VehiclePara.Iz;
Ms =1430;
Np = MPCParameters.Np;
Ns = MPCParameters.Ns;
Nx = MPCParameters.Nx;
Ny = MPCParameters.Ny;
Ts = MPCParameters.Ts;
Tsl = MPCParameters.Tsl; 
dy  = -lr;
% Np      = 25; % 0.1
% Ns      = 10; % 0.05;
% Nx      = 8;
% Ny      = 5;
% Ts      = 0.01;
% Tsl     = 0.5;

%% 
K_phi       = 145330;  % roll stiffness coefficient
D_phi       = 4500;    % roll damping coefficient
Theta_1     = Calpha_f + Calpha_r;  
Theta_2     = lf*Calpha_f - lr*Calpha_r;
Theta_3     = lf*lf*Calpha_f + lr*lr*Calpha_r;
Theta_4     =  Ix*M - Ms*hCG ;
Theta_5     =  Ms*hCG*g - K_phi;
Theta_6     =  Ms*hCG;

Ixz     = 0; %500
Mint = [  M         0       -M*hCG          0;
          0         Iz      -Ixz            0;
          -M*hCG    -Ixz    Ix+M*hCG*hCG    0
          0         0       0               1];

Nint = [-Theta_1/Vel    M*Vel-Theta_2/Vel   0           M*g;  %0-->M*g
        -Theta_2/Vel    -Theta_3/Vel        0           0;
        0               -M*hCG*Vel          D_phi       K_phi-M*g*hCG;
        0               0                   -1          0];

F1int = [-Calpha_f; -lf*Calpha_f;  0;  0];
F2int = [   0        0; 
            0        0;
            K_phi    0; %-M*g*hCG
            0        0];

Ac_11     = -Mint\Nint; % 4*4
B1cn_11   = Mint\F1int; % 4*1
B2cn_11   = Mint\F2int; % 4*2

Ac_12     = [0 0; 0 0; 0 0; 0 0];
Ac_21     = [1  0    0    0;
             0  1    0    0];
Ac_22     = [0   Vel; 0   0]; 
 
B1cn_12   = [ 0; 	0]; 
B2cn_12   = [ 0, 	0; 
              0, 	-Vel];

Acn     = [Ac_11 Ac_12; Ac_21 Ac_22];
B1cn    = [B1cn_11; B1cn_12];
B2cn    = [B2cn_11; B2cn_12];

Theta_1     = Calpha_f + Calpha_r;  
Theta_2     = lf*Calpha_f - lr*Calpha_r;
Theta_3     = lf*lf*Calpha_f + lr*lr*Calpha_r;
Theta_4     =  Ix*M - Ms*hCG ;
Theta_5     =  Ms*hCG*g - K_phi;
Theta_6     =  Ms*hCG;

Mint = [  M*Vel      0       -M*hCG         0;
          0          Iz      -Ixz           0;
          M*Vel*hCG  -Ixz    Ix+Theta_6*Theta_6  0;
          0          0       0               1];

Nint = [Theta_1    -M*Vel + Theta_2/Vel     0     M*g;  %0-->M*g
        Theta_2    Theta_3/Vel              0       0;
        0            -M*hCG*Vel             Theta_5   D_phi;
        0               0                   -1      0 ];

F1int = [-Calpha_f; -lf*Calpha_f;  0;  0];
F2int = [   0        0; 
            0        0;
            -K_phi   0; %-M*g*hCG
            0        0];

Ac_11     = Mint\Nint; % 4*4
B1cn_11   = Mint\F1int; % 4*1
B2cn_11   = Mint\F2int; % 4*2

A_22 =  2*(Calpha_f + Calpha_r)/ Vel/ M;
A_23 = -2*(Calpha_f +Calpha_r)/ M;
A_24 =  2*(lf * Calpha_f - lr * Calpha_r)/ Vel/ M;
A_42 = (lf * Calpha_f - lr * Calpha_r)/ Vel/ Iz;
A_43 = -(lf * Calpha_f - lr * Calpha_r)/ Iz;
A_44 = (lf^2 * Calpha_f + lr^2 * Calpha_r)/ Vel/ Iz;
% A_e = [0 1 0 0;
%         0 A_22 A_23 A_24;
%         0 0 0 1;
%         0 A_42 A_43 A_44];
% B_e = [0; -Calpha_f/M; 0; -lf*Calpha_f/Iz];
ac_12   = zeros(4,4);
ac_21   = zeros(4,4);       
ac_22   =  [0 1 0 0;
            0 A_22 A_23 A_24;
            0 0 0 1;
            0 A_42 A_43 A_44];
acn      = [Ac_11 ac_12; ac_21 ac_22];

b1cn_13 = [0; -Calpha_f/M; 0; -lf*Calpha_f/Iz];
b2cn_13 = [0; A_24*Vel - Vel; 0; A_44 * Vel];

b1cn    = [B1cn_11; b1cn_13];
b2cn    = [B2cn_11; zeros(4,1), b2cn_13];

N_x = size(acn,2);
N_u = size(b1cn,2);
N_w = size(b2cn,2);
Mc = [acn, b1cn, b2cn;
      zeros(N_u + N_w, N_x + N_u + N_w)];
Md = expm(Ts * Mc);
acd1 = Md(1:N_x, 1:N_x);
bcd1 = Md(1:N_x, N_x + N_u); 
bcd2 = Md(1:N_x, N_x + N_u + 1 : N_x + N_u + N_w);
% Linear bicycle model
% State: [vy, r, ey, epsi]
% Input: [delta_f];
% reference: [Kappa]
Ab = zeros(2,2);

Ab(1,1) =  (Calpha_f + Calpha_r)/ Vel/ M;
Ab(1,2) =  (lf * Calpha_f - lr * Calpha_r)/ Vel/ M - Vel;
Ab(2,1) =  (lf * Calpha_f - lr * Calpha_r)/ Vel/ Iz;
Ab(2,2) =  (lf^2 * Calpha_f + lr^2 * Calpha_r)/ Vel/ Iz;

Bb = zeros(2,1);
Bb(1,1) = -Calpha_f/ M; 
Bb(2,1) = -lf* Calpha_f/ Iz;

Au = [Ab, zeros(2,2);
     1, 0, 0, Vel;
     0, 1,  0, 0;];

Bu = [Bb; 0; 0];
Br = [0; zeros(2,1); 0];

N_x = size(Au,2);
N_u = size(Bu,2);
N_w = size(Br,2);
Mc = [Au, Bu, Br;
      zeros(N_u + N_w, N_x + N_u + N_w)];
Md = expm(Ts * Mc);
Ad1 = Md(1:N_x, 1:N_x);
Bd1 = Md(1:N_x, N_x + N_u); 
Bd2 = Md(1:N_x, N_x + N_u + N_w);

%% lqr for error-based model
%% states as [ed_dot, ed, epsi_dot, epsi]
%% input as delta_f
%% states obtained by the vehicle position and the planned position
Calpha_f = VehiclePara.CafHat;
Calpha_r = VehiclePara.CarHat;
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

B_r = [0; (lf * Calpha_f - lr * Calpha_r)/ M - Vel^2; 0; (lf^2 * Calpha_f + lr^2 * Calpha_r)/ Iz];


%% ZOH/FOH
    As = eye(Nx) + Ts * Acn;
    Bs1 = Ts * B1cn; 
    Bs2 = Ts * B2cn;

    As4 = eye(4) + Ts * Au;
    Bs14 = Ts * Bu; 
    Bs24 = Ts * Br;

    AcnTsl = Tsl*Acn;
    Al = eye(Nx) + AcnTsl + 0.5*AcnTsl*AcnTsl  + (AcnTsl*AcnTsl*AcnTsl)/6; %; % 

    Baug    = [B1cn, B2cn];
    BaugTsl = Baug*Tsl;
    Gail1   = BaugTsl + 0.5*(AcnTsl * BaugTsl)+ (AcnTsl * AcnTsl * BaugTsl)/6;%; % 
    Gail2   = 0.5*BaugTsl + (AcnTsl * BaugTsl)/6;%;%
    
    B1l     = Gail1 - Gail2; 
    B2l     = Gail2;
    Bl11    = B1l(:,1);
    Bl12    = B1l(:,2:3);
    Bl21    = B2l(:,1);
    Bl22    = B2l(:,2:3);

    StateSpaceModel.AC      = Ac_11;
    StateSpaceModel.BC1     = B1cn_11;
    StateSpaceModel.BC2     = B2cn_11;
    StateSpaceModel.Acn     = Acn;
    StateSpaceModel.B1cn    = B1cn; 
    StateSpaceModel.B2cn    = B2cn; 

    StateSpaceModel.acn     = acn;
    StateSpaceModel.b1cn    = b1cn; 
    StateSpaceModel.b2cn    = b2cn;  
    
    StateSpaceModel.acd1     = acd1;
    StateSpaceModel.bcd1     = bcd1; 
    StateSpaceModel.bcd2     = bcd2; 
    
    % intitial portion of the prediction horizon
    StateSpaceModel.As      = As;
    StateSpaceModel.Bs1     = Bs1; 
    StateSpaceModel.Bs2     = Bs2;
    
    % latter portion of the prediction horizon
    StateSpaceModel.Al      = Al;    % state matrix at timestep k
    StateSpaceModel.Bl11    = Bl11;  % inputs matrix at timestep k
    StateSpaceModel.Bl12    = Bl12;  % disturbance matrix at timestep k
    StateSpaceModel.Bl21    = Bl21;  % inputs matrix at timestep k+1
    StateSpaceModel.Bl22    = Bl22;  % disturbance matrix at timestep k+1

    % new state space model

    StateSpaceModel.As4     = As4;
    StateSpaceModel.Bs14    = Bs14;
    StateSpaceModel.Bs24    = Bs24;

    StateSpaceModel.Anom    = Ab;
    StateSpaceModel.Bnom    = Bb;

    StateSpaceModel.Au      = Au;
    StateSpaceModel.Bu      = Bu;
    StateSpaceModel.Br      = Br;

    StateSpaceModel.A_e     = A_e;
    StateSpaceModel.B_e     = B_e;

    StateSpaceModel.Ad1     = Ad1;
    StateSpaceModel.Bd1     = Bd1;
    StateSpaceModel.Bd2     = Bd2;
end % end of func_SpatialDynamicalModel

