function [sys,x0,str,ts] = II1(t,x,u,flag)
%***************************************************************%
    switch flag,
        case 0 % Initialization %
            [sys,x0,str,ts] = mdlInitializeSizes; % Initialization
        case 2 % Update %
            sys = mdlUpdates(t,x,u); % Update discrete states
        case 3 % Outputs %
            sys = mdlOutputs(t,x,u); % Calculate outputs
        case {1,4,9} % Unused flags
            sys = [];            
        otherwise % Unexpected flags %
            error(['unhandled flag = ',num2str(flag)]); % Error handling
    end %  end of switch    
%  End sfuntmpl

function [sys,x0,str,ts] = mdlInitializeSizes
%==============================================================
% Initialization, flag = 0��mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%==============================================================
sizes = simsizes;
sizes.NumContStates  = 0; 
sizes.NumDiscStates  = 6;  
sizes.NumOutputs     = 43; 
sizes.NumInputs      = 38; 
sizes.DirFeedthrough = 1;  
sizes.NumSampleTimes = 1; 
sys = simsizes(sizes);  
x0 = zeros(sizes.NumDiscStates,1);

str = [];          

ts  = [0.05 0];       % ts=[period, offset].sample time=0.05,50ms 

%--Global parameters and initialization
[y, e] = func_RLSFilter_Calpha_f('initial', 0.95, 10, 10);
[y, e] = func_RLSFilter_Calpha_r('initial', 0.95, 10, 10);

    global InitialGapflag; 
    InitialGapflag = 0; % the first few inputs don't count. Gap it.
    
    global VehiclePara; % for SUV
    VehiclePara.m   = 1460;   %m,Kg; Sprung mass = 1370
    VehiclePara.g   = 9.81;
    VehiclePara.hCG = 0.65;%m
    VehiclePara.Lf  = 1.05;  % 1.05
    VehiclePara.Lr  = 1.55;  % 1.55
    VehiclePara.L   = 2.6;  %VehiclePara.Lf + VehiclePara.Lr;
    VehiclePara.Tr  = 1.565;  %c,or 1.57. 
    VehiclePara.mu  = 0.85; % 0.55;
    VehiclePara.Iz  = 2059.2;     
    VehiclePara.Ix  = 700.7;   
    VehiclePara.Radius = 0.387;  
    
    global MPCParameters; 
    MPCParameters.Np  = 5;% predictive horizon Assume Np=Nc
    MPCParameters.Ns  = 5; %  Tsplit
    MPCParameters.Ts  = 0.05; % the sample time of near term
    MPCParameters.Tsl = 0.05; % the sample time of long term       
    MPCParameters.Nx  = 6; %the number of state variables
    MPCParameters.Ny  = 2; %the number of output variables      
    MPCParameters.Nu  = 1; %the number of control inputs
    
    global CostWeights; 
    CostWeights.Wephi   = 10; %state vector =[beta,yawrate,e_phi,s,e_y]
    CostWeights.Wey     = 100;
    CostWeights.Weyd    = 5;
    CostWeights.Wepsid  = 5;
    CostWeights.Ddeltaf = 1;
    CostWeights.deltaf  = 1;
    CostWeights.Wshar   = 500;
    CostWeights.Wshr    = 500;
    CostWeights.LTR     = 500;
    CostWeights.ZMP     = 500;

    global Constraints;  
    Constraints.dumax   = 0.1/MPCParameters.Ts;     % Units: rad,0.08rad/s--4.6deg/s  
    Constraints.umax    = 0.4;      % Units: rad appro.23deg
    Constraints.dumax_L  = 0.08;
    Constraints.arlim   = 6*pi/180; % alpha_lim=6deg~ 0.1047rad
    Constraints.rmax    = 1.0; % rr_max = 1rad/s    
    Constraints.LRT     = 0.9;
    Constraints.zmp     = 0.7;

    ar_slackMax         = 6*pi/180; % rad
    rmax_slackMax       = 1.0;
    Constraints.Sshmax  = [ar_slackMax; rmax_slackMax];
    
    Constraints.DPhimax = 90*pi/180;  % 
    Constraints.Dymax   = 1.7; % 3.0;    cross-track-error max 3m

    global WayPoints_IndexPre;
    WayPoints_IndexPre = 1;
    
    global Reftraj;
    % Reftraj = load('WayPoints_Alt3fromFHWA_Samples.mat');
    Reftraj = load('PathPoints_DoubleLane_new.mat');    
    
    global FeedbackCorrects;
    FeedbackCorrects.StatePred = zeros(6,1);
    FeedbackCorrects.Ctrlopt   = 0;
    
    Vel = 60/3.6;
    VehiclePara.CarHat = -77000;
    VehiclePara.CafHat = -82000;
    global StateSpaceModel
    [StateSpaceModel] = func_RigidbodyModel_FOH_Matrix_ROLL_YALMIP(VehiclePara, MPCParameters, Vel, VehiclePara.CafHat, VehiclePara.CarHat);
    
    Np          = MPCParameters.Np;
    Eymax       = zeros(Np,1);
    Eymin       = zeros(Np,1);     
    LM_right    = -5;
    LM_middle   = 0;
    Yroad_L     = -2.5;
    for i =1:1:Np  % ע��ey�Ǵ����ŵ�, Np = 25
        Eymax(i,1) = (LM_middle - Yroad_L);
        Eymin(i,1) = (LM_right - Yroad_L);             
    end

    [Envelope] = func_SafedrivingEnvelope_SL(VehiclePara, MPCParameters, Constraints, StateSpaceModel, Vel, VehiclePara.CarHat, Eymax, Eymin); 
    [Sl, Ql, Rdun,Rdul, Wshl, dun, dul] = func_CostWeightingRegulation_QuadSlacks_Yalmip(MPCParameters, CostWeights, Constraints);
   
    CostWeights.Ql = Ql;    % Weights on states
    CostWeights.Sl = Sl;    % Weight  on fwa
    CostWeights.Rdun = Rdun; % Weight on dfwa
    CostWeights.Rdul = Rdul;
    CostWeights.Wshl = Wshl; % Weight on slack variables
    CostWeights.dun  = dun;
    CostWeights.Q_err  = dul;
    global RMPC
    global solution
    global fwa_sol
    global oldu 

    oldu = zeros(1,MPCParameters.Np);

    RMPC = Tube_MPC_ROLL3(VehiclePara, MPCParameters, Constraints, CostWeights, StateSpaceModel, Envelope, Vel);

    Q_e = diag([0.5 0.1 1 1]);
    R_e   = 1;
    A_e = StateSpaceModel.A_e;
    B_e = StateSpaceModel.B_e;
    K = lqr(A_e, B_e, Q_e, R_e);
    StateSpaceModel.K = K;
    solution = 0;
    fwa_sol  = 0;
%  End of mdlInitializeSizes

function sys = mdlUpdates(t,x,u)
%==============================================================
% Update the discrete states, flag = 2�� mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%==============================================================
    sys = x;    
% End of mdlUpdate.

function sys = mdlOutputs(t,x,u)
%==============================================================
% Calculate outputs, flag = 3�� mdlOutputs
% Return the block outputs. 
%==============================================================

%***********Step (1). Parameters Initialization ***************************************%

global InitialGapflag;
global VehiclePara;
global MPCParameters; 
global CostWeights;     
global Constraints;
global WayPoints_IndexPre;
global Reftraj;
global FeedbackCorrects;
global RMPC
global oldu
global solution
global fwa_sol
global X_nominal
global StateSpaceModel

Ctrl_SteerSW    = 0;
t_Elapsed       = 0;
PosX            = 0;
PosY            = 0;
PosPsi          = 0;
PSI             = 0;
Vel             = 0;
e_psi           = 0;
e_y             = 0;
fwa_opt         = 0;
Shenvelop_hat   = [0; 0];
r_ssmax         = 0;
YZPM            = 0; 
y_zmp           = 0; 
LTR             = 0; 
LTR_T           = 0;
ENV             = [0;0];
Vy              = 0;
alphar          = 0;
Roll            = 0;
Roll_Shad       = 0;
Roll_BaknR      = 0;
Station         = 0;
yawrate         = 0;
Beta            = 0;
Nx              = 8;
Np              = MPCParameters.Np;
x0              = zeros(Nx,1);
x_nom1          = zeros(Nx,1);
x_nom2          = zeros(Nx,1);
x_nom3          = zeros(Nx,1);
dist            = zeros(Nx, Np);
Dist            = zeros(Np *2, 1);
epsi_mean       = 0;
fwa_du_1        = 0;
fwa_du_sol      = zeros(1,Np);
fwa_current     = 0;
fwa_e           = 0;
err             = zeros(4,1);
forwa           = 0;
U               = zeros(2,Np);
cost            = 0;
CafHat      = 0;
CarHat      = 0;
Fyf         = 0;
Fyr         = 0;
Arfa_f      = 0;
Arfa_r      = 0;
Ay          = 0;
tire_parameter = zeros(6,1);
if InitialGapflag < 2 %  get rid of the first two inputs,  because no data from CarSim
    InitialGapflag = InitialGapflag + 1;
else % start control
    InitialGapflag = InitialGapflag + 1;
%***********Step (2). State estimation and Location **********************% 
    t_Start = tic; % ��ʼ��ʱ  
    %-----Update State Estimation of measured Vehicle Configuration--------%
    [VehStateMeasured, ParaHAT] = func_StateEstimation(u);   
    PosX        = VehStateMeasured.X;
    PosY        = VehStateMeasured.Y;
    PosPsi      = VehStateMeasured.phi;    
    Vel         = VehStateMeasured.x_dot; 
    Vy          = VehStateMeasured.y_dot; 
    yawrate     = VehStateMeasured.phi_dot; % rad/s
    Ax          = VehStateMeasured.Ax; % x_dot
    Ay          = VehStateMeasured.Ay; % y_dot

    delta_f     = VehStateMeasured.delta_f;% deg-->rad    
%     Steer_SW    = VehStateMeasured.Steer_SW;
    fwa         = VehStateMeasured.fwa;
    Beta        = VehStateMeasured.beta;%rad
    Roll_Shad   = VehStateMeasured.Roll_Shad;%rad
    Station     = VehStateMeasured.Station;
    Roll        = ParaHAT.Roll;
    Rollrate    = ParaHAT.Rollrate;
%     Ay_CG       = ParaHAT.Ay_CG;    
%     Ay_Bf_SM    = ParaHAT.Ay_Bf_SM;    
    Fyf         = ParaHAT.Fyf;
    Fyr         = ParaHAT.Fyr;   
    alphaf      = ParaHAT.alphaf;
    alphar      = ParaHAT.alphar;
    
    %-----Estimate Cornering stiffness -------------------%  
%     for front tire
    Arfa_f = (Beta + yawrate*VehiclePara.Lf/Vel - fwa_sol);
    [yf, Calpha_f1] = func_RLSFilter_Calpha_f(Arfa_f, Fyf);
    CafHat = sum(Calpha_f1);
    % if CafHat > -30000
    %     CafHat = -70000;
    % end
    %for rear tire 
    Arfa_r = (Beta - yawrate*VehiclePara.Lr/Vel);
    [yr, Calpha_r1] = func_RLSFilter_Calpha_r(Arfa_r, Fyr);
    CarHat = sum(Calpha_r1);
    % if CarHat > -30000
    %     CarHat = -92000;
    % end    
    % %-----OR use constant Cornering stiffness -------------------%  
    % CafHat = -45000;
    % CarHat = -45000;
    
    %********Step(3): Given reference trajectory, update vehicle state and bounds *******************% 
    [WPIndex, ~, ~, Uaug, ~, PrjP, ~] = func_RefTraj_LocalPlanning_TwoTimeScales_Spatial_Integrated( MPCParameters,... 
                            VehiclePara,... 
                            WayPoints_IndexPre,... 
                            Reftraj.DL_PP,... 
                            VehStateMeasured); %                                     
    if ( WPIndex <= 0)
       fprintf('Error: WPIndex <= 0 \n');    %����
    else
        Xm = [Vy; yawrate; Rollrate; Roll; PrjP.ey; PrjP.epsi];        
        WayPoints_IndexPre = WPIndex;        
    end

    %****Step(4):  MPC formulation;********************%
    Station =  PrjP.yr;  
    Np          = MPCParameters.Np;
    Eymax       = zeros(Np,1);
    Eymin       = zeros(Np,1);     
    LM_right    = -5;
    LM_middle   = 0;
    Yroad_L     = -2.5;
    for i =1:1:Np 
        Eymax(i,1) = (LM_middle - Yroad_L);
        Eymin(i,1) = (LM_right - Yroad_L);             
    end
    [Envelope] = func_SafedrivingEnvelope_SL(VehiclePara, MPCParameters, Constraints, StateSpaceModel, Vel, CarHat, Eymax, Eymin);

    U = reshape(cell2mat(Uaug), 2,  MPCParameters.Np);
%% obtain the steer_wheel_angle from LQR (feedback + feedforward)
    Kappar = U(2,1);
    K = StateSpaceModel.K;
    err1  = PrjP.ey;
    err2  = Vy * cos(PrjP.epsi) + Vel * sin(PrjP.epsi);
    err3  = PrjP.epsi;
    ss_dot = -Vy * sin(PrjP.epsi_m) + Vel * cos(PrjP.epsi);
    err4  = yawrate - Kappar * Vel;
    err   = [err1; err2; err3; err4];
    fwa_e = -K * err;
%% obtain the SWA from Stanley method
        x0 = [Beta; yawrate; Rollrate; Roll_Shad; err];
        U = reshape(cell2mat(Uaug), 2,  MPCParameters.Np);
        fwa_current = oldu(solution+1);
        inputs = {x0, U, fwa_current};
        
        [solutions, diagnostics] = RMPC{inputs};

        if diagnostics
        solution    = solution + 1;
        fwa_sol     = oldu(solution);
        else
        % nominal states in the controller
        x_nom       = solutions{1};
        x_nom1      = x_nom(:,2);
       fwa_solution = solutions{2};
       oldu         = fwa_solution;
       fwa_sol      = fwa_solution(1);
       solution     = 1;             
       fwa_du_sol   = solutions{3};
        cost         = solutions{4};
        end

    % end
    %====================================================================%
    Ctrl_SteerSW = 19 * fwa_sol*180/pi; % in deg.    
      
    t_Elapsed = toc( t_Start ); %computation time
    
    %-----------------------------------------%
    e_y            = PrjP.ey;
    e_psi          = PrjP.epsi; 
    PSI            = PrjP.Psi_10;

    Bs2            = StateSpaceModel.bcd2;
    dist           = Bs2 * U;
    Dist           = [dist(6,:)';dist(8,:)'];

    tire_parameter = [CafHat; CarHat; Fyf; Fyr; Arfa_f; Arfa_r];

    x_r = [Vy; yawrate; Rollrate; Roll_Shad];

    Hsh  = [1/Vel   -VehiclePara.Lr/Vel    0   0    ;
            0         1                    0   0    ];   
    Psh  = [0  0;  VehiclePara.g/Vel   0];
    ENV  = Hsh * x_r + Psh * U(:,1);
    Fzl = ParaHAT.Fz_l1 + ParaHAT.Fz_l2;
    Fzr = ParaHAT.Fz_r1 + ParaHAT.Fz_r2;
    LTR = (Fzr - Fzl)./(Fzr + Fzl);
    K_phi       = 145330;  % roll stiffness coefficient
    D_phi       = 4500;    % roll damping coefficient
    H_LTR    = 2 / (VehiclePara.m*VehiclePara.g * VehiclePara.Tr) * [ 0, 0, D_phi, K_phi];
    LTR_T = H_LTR * x_r;
end % end of if Initialflag < 2 % 
sys = [Ctrl_SteerSW; t_Elapsed; PosX; PosY; PosPsi; Station; PSI; e_psi; e_y; cost; LTR_T; LTR; Vy; Beta; yawrate; fwa_sol; fwa_e; err; U(2,1:end)'; solution;Dist; tire_parameter]; %

%***************************************************************%
% **** State estimation
%***************************************************************%
function [VehStatemeasured, HATParameter] = func_StateEstimation(ModelInput)
%***************************************************************%
% we should do state estimation, but for simplicity we deem that the
% measurements are accurate
% Update the state vector according to the input of the S function,
%           usually do State Estimation from measured Vehicle Configuration
%***************************************************************%  
    g = 9.81;
    VehStatemeasured.X       = round(100*ModelInput(1))/100;
    VehStatemeasured.Y       = round(100*ModelInput(2))/100;%  
    VehStatemeasured.phi     = (round(10*ModelInput(3))/10)*pi/180; 
    VehStatemeasured.x_dot   = ModelInput(4)/3.6; %Unit:km/h-->m/s 
    VehStatemeasured.y_dot   = ModelInput(5)/3.6; %Unit:km/h-->m/s  
    VehStatemeasured.phi_dot = (round(10*ModelInput(6))/10)*pi/180; %Unit��deg/s-->rad/s
    VehStatemeasured.beta    = (round(10*ModelInput(7))/10)*pi/180;% side slip, Unit:deg-->rad
    VehStatemeasured.delta_f = (round(10*0.5*(ModelInput(8)+ ModelInput(9)))/10)*pi/180; % deg-->rad
    VehStatemeasured.delta_f = 0.5*(ModelInput(8)+ ModelInput(9))*pi/180; % deg-->rad

    VehStatemeasured.fwa     = VehStatemeasured.delta_f * pi/180;  % deg-->rad
    VehStatemeasured.Steer_SW= ModelInput(10); %deg
    VehStatemeasured.Ax      = g*ModelInput(11);%
    VehStatemeasured.Ay      = g*ModelInput(12);%
    VehStatemeasured.yawrate_dot = ModelInput(13); %rad/s^2
    % VehStatemeasured.SlipAngle = ModelInput(39);
    % Here I don't explore the state estimation process, and deem the
    % measured values are accurate!!! 
    HATParameter.alpha_l1   = (round(10*ModelInput(14))/10)*pi/180; % deg-->rad
    HATParameter.alpha_l2   = (round(10*ModelInput(15))/10)*pi/180; % deg-->rad 
    HATParameter.alpha_r1   = (round(10*ModelInput(16))/10)*pi/180; % deg-->rad
    HATParameter.alpha_r2   = (round(10*ModelInput(17))/10)*pi/180; % deg-->rad  
    HATParameter.alphaf     = (round(10*0.5 * (ModelInput(14)+ ModelInput(16)))/10)*pi/180; % deg-->rad 
    HATParameter.alphar     = (round(10*0.5 * (ModelInput(15)+ ModelInput(17)))/10)*pi/180; % deg-->rad
    
    HATParameter.Fz_l1      = round(10*ModelInput(18))/10; % N 
    HATParameter.Fz_l2      = round(10*ModelInput(19))/10; % N 
    HATParameter.Fz_r1      = round(10*ModelInput(20))/10; % N 
    HATParameter.Fz_r2      = round(10*ModelInput(21))/10; % N 
    
    HATParameter.Fy_l1      = round(10*ModelInput(22))/10; % N 
    HATParameter.Fy_l2      = round(10*ModelInput(23))/10; % N 
    HATParameter.Fy_r1      = round(10*ModelInput(24))/10; % N 
    HATParameter.Fy_r2      = round(10*ModelInput(25))/10; % N 
    HATParameter.Fyf        = HATParameter.Fy_l1 + HATParameter.Fy_r1;
    HATParameter.Fyr        = HATParameter.Fy_l2 + HATParameter.Fy_r2;
    
    HATParameter.Fx_L1      = ModelInput(26);
    HATParameter.Fx_L2      = ModelInput(27);
    HATParameter.Fx_R1      = ModelInput(28);
    HATParameter.Fx_R2      = ModelInput(29);
    
%     HATParameter.GearStat    = ModelInput(30);
    VehStatemeasured.Roll_Shad   = ModelInput(30)*pi/180;% deg-->rad 
    HATParameter.Roll        = ModelInput(31)*pi/180;% deg-->rad 
    HATParameter.Rollrate    = ModelInput(32)*pi/180;% deg/s-->rad/s
    HATParameter.Roll_accel  = ModelInput(33); % rad/s^2
    HATParameter.Z0          = ModelInput(34); %m
    VehStatemeasured.Station     = ModelInput(35); %m
    HATParameter.Zcg_TM      = ModelInput(35); %m
    HATParameter.Zcg_SM      = ModelInput(36); %m
    HATParameter.Ay_CG       = ModelInput(37)*g; %m/s^2
    HATParameter.Ay_Bf_SM    = ModelInput(38)*g; %m/s^2
    
% end % end of func_StateEstimation
