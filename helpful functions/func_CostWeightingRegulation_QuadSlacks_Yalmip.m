function [Sl, Ql, Rdun,Rdul, Wshl, dun, Q_err] = func_CostWeightingRegulation_QuadSlacks_Yalmip(MPCParameters, CostWeights, Constraints)

%% ������ʼ��
    Ts      = MPCParameters.Ts;
    Tsl     = MPCParameters.Tsl;
    
    Qephi   = CostWeights.Wephi;
    Qey     = CostWeights.Wey;
    Qephid  = CostWeights.Wepsid;
    Qeyd    = CostWeights.Weyd;
    RDdeltaf= CostWeights.Ddeltaf;   
    Sdeltaf = CostWeights.deltaf;    
    Wshar   = CostWeights.Wshar;
    Wshr    = CostWeights.Wshr;

    DPhimax = Constraints.DPhimax;  %  0.15 rad ==> 8.5deg
    Dymax   = Constraints.Dymax;
    armax   = Constraints.arlim; % 0.104rad=6deg; %0.15rad=8deg
    rmax    = Constraints.rmax;    
    
    umax    = Constraints.umax;
    dumax   = Constraints.dumax;
    dun     = dumax * Ts;
    dul     = dumax * Tsl;
    
    %% Ȩ�����ӹ�һ��
    Qephi_DPhimax2  = Qephi/(DPhimax*DPhimax);
    Qey_Dymax2      = Qey/(Dymax*Dymax);    
    Ql              = diag([0.5, 0.5, Qey_Dymax2, Qeyd, Qephi_DPhimax2, Qephid]);
    Q_err           = diag([Qey_Dymax2, Qeyd, Qephi_DPhimax2, Qephid]);
    Wshr_rmax2      = Wshr/(rmax*rmax);    
    Wshar_armax2    = Wshar/(armax*armax);
    Wshl            = diag([Wshar_armax2, Wshr_rmax2]);
   
    dumax_ts2   = (dun * dun);%
    dumax_tl2   = (dul * dul);% 
    Rdun        = RDdeltaf/dumax_ts2;
    Rdul        = diag([Qeyd, Qephid, Qey_Dymax2, Qephi_DPhimax2]);
    
    Sl = Sdeltaf/(umax * umax);
    
end % end of func_CostWeightingRegulation