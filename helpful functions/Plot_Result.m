% close all
lw = 2;
Num = length(u.signals.values(:,2));
output          = u.signals.values;
SteerWheel      = output(:,1);
time_elapsed    = output(:,2);
PosX            = output(:,3);
PosY            = output(:,4);
PosPsi          = output(:,5);
Station         = output(:,6);
PSI_reference   = output(:,7);
Error_psi       = output(:,8);
Error_ey        = output(:,9);
Objective       = output(:,10);
LTR_t           = output(:,11:12);
Vy              = output(:,13);
Sideslip        = output(:,14);
Yawrate         = output(:,15);
Actual_input    = output(:,16);
Nominal_input   = output(:,17);
Error           = output(:,18:21);
Envelope        = output(:,22:23);
Nominal_state   = output(:,24:31);
Actual_state    = output(:,32:39);
De_Input        = output(:,40:44);
Curvature       = output(:,45:48);
Sol_flag        = output(:,49);
Tire_state      = output(:,50:55);

%% Ey_dot---Epsi_dot
Ey_dot_nomi     = Nominal_state(:,6);
Epsi_dot_nomi   = Nominal_state(:,8);
X1(1,:)         = Ey_dot_nomi;
X1(2,:)         = Epsi_dot_nomi;
Ey_dot_actu     = Actual_state(:,6);
Epsi_dot_actu   = Actual_state(:,8);
X2(1,:)         = Ey_dot_actu;
X2(2,:)         = Epsi_dot_actu;
De_X            = X2 - X1;
%% plot the trajectory of the state for both actual and nominal
Num_ini = 606;
figure(1);
plot(X1(1,Num_ini),X1(2,Num_ini),'--k+');
plot(X2(1,Num_ini),X2(2,Num_ini),'--b*');
%axis([-20,20,-10,10]);
hold on;
%%
Z=[z1(1)+X1(1,Num_ini),z2(1)+X1(1,Num_ini),-z1(1)+X1(1,Num_ini),-z2(1)+X1(1,Num_ini),z1(1)+X1(1,Num_ini);
    z1(2)+X1(2,Num_ini),z2(2)+X1(2,Num_ini),-z1(2)+X1(2,Num_ini),-z2(2)+X1(2,Num_ini),z1(2)+X1(2,Num_ini)];
plot(Z(1,:),Z(2,:), 'k');
h0=fill(Z(1,:),Z(2,:), [0.7,0.8,1]);
set(h0,'edgealpha',0,'facealpha',0.3) 
%%
for j= Num_ini: Num_ini + 10
Z=[z1(1)+X1(1,j),z2(1)+X1(1,j),-z1(1)+X1(1,j),-z2(1)+X1(1,j),z1(1)+X1(1,j);
    z1(2)+X1(2,j),z2(2)+X1(2,j),-z1(2)+X1(2,j),-z2(2)+X1(2,j),z1(2)+X1(2,j)];
% Z=[z1(1)+X1(1,j),z2(1)+X1(1,j),z3(1)+X1(1,j),-z1(1)+X1(1,j),-z2(1)+X1(1,j),-z3(1)+X1(1,j),z1(1)+X1(1,j);
%     z1(2)+X1(2,j),z2(2)+X1(2,j),z3(2)+X1(1,j),-z1(2)+X1(2,j),-z2(2)+X1(2,j),-z3(2)+X1(1,j),z1(2)+X1(2,j)];
plot(Z(1,:),Z(2,:),'k');
h0=fill(Z(1,:),Z(2,:), [0.7,0.8,1]);
set(h0,'edgealpha',0,'facealpha',0.3) 
end
%%
plot(X1(1,Num_ini:Num_ini+20) ,X1(2,Num_ini:Num_ini+20),'--ro','MarkerFaceColor','r');
plot(X2(1,Num_ini:Num_ini+20) ,X2(2,Num_ini:Num_ini+20),'--bs','MarkerFaceColor','b');
hold on
%%
i = Num_ini + 20;
plotregion([P;-P],[-b;-b],[],[],'r');hold on
plot(X1(1,i),X1(2,i),'--ro','MarkerFaceColor','r');
plot(X2(1,i),X2(2,i),'--bs','MarkerFaceColor','b');
grid on
%%
%plot 3D graph for states trajectory
figure(2);
%plot the trajectory of the state for both actual and nominal
k=[];
start=Num_ini;
for j=start:i
Z=[-z1(1)+X1(1,j),-z2(1)+X1(1,j),z1(1)+X1(1,j),z2(1)+X1(1,j),-z1(1)+X1(1,j);
   -z1(2)+X1(2,j),-z2(2)+X1(2,j),z1(2)+X1(2,j),z2(2)+X1(2,j),-z1(2)+X1(2,j)];
plot3([j*ones(1,5)],Z(1,:),Z(2,:),'r')
k=[k,j];
hold on;
end
grid on;
hold on;
plot3(k,X1(1,start:i),X1(2,start:i),'-k+','linewidth',1);
plot3(k,X2(1,start:i),X2(2,start:i),'-b*','linewidth',1);
title('Tube-Based MPC');
ylabel('$x_1$','interpreter','latex')
zlabel('$x_2$','interpreter','latex')
xlabel('Number of Interation','interpreter','latex')
%%

figure(3);
%plot the trajectory of the state for both actual and nominal 
k=[]; i= 20; X1 = X1(:,Num_ini +1 :i+Num_ini);X2 = X2(:,Num_ini +1 :i+Num_ini);
for j=1:i
Z=[-z1(1)+X1(1,j),-z2(1)+X1(1,j),z1(1)+X1(1,j),z2(1)+X1(1,j),-z1(1)+X1(1,j);
    -z1(2)+X1(2,j),-z2(2)+X1(2,j),z1(2)+X1(2,j),z2(2)+X1(2,j),-z1(2)+X1(2,j)];

plot3([j*ones(1,5)],Z(1,:),Z(2,:),'r')
fill3(j*ones(1,5),Z(1,:),Z(2,:), [1,1,0]);
k=[k,j];
hold on;
end
grid on;

plot3(1,X1(1,1),X1(2,1),'k+');
plot3(1,X2(1,1),X2(2,1),'r*');
hold on;

plot3(k,X1(1,:),X1(2,:),'-k+','linewidth',1);
plot3(k,X2(1,:),X2(2,:),'-r*','linewidth',1);
title('Tube-Based MPC');
ylabel('$x_1$','interpreter','latex')
zlabel('$x_2$','interpreter','latex')
xlabel('Number of Interation','interpreter','latex')

%%
% LTR = u.signals.values(:,9);
% Yzmp_SH = u.signals.values(:,10);
% Yzmp = u.signals.values(:,11);
% Yzmp1 = u.signals.values(:,12);

% figure
% plot(1:Num, LTR','r',1:Num, Yzmp_SH,'b',1:Num, Yzmp,'g',1:Num, Yzmp1,'k', 'Linewidth',lw);
% grid on
% legend('LTR','Yzmp_{SH}','Yzmp','Yzmp1');

%%

% Yzmp_SH = u.signals.values(:,10);
% Yzmp = u.signals.values(:,11);
% LTR = u.signals.values(:,12);

% figure
% plot(1:Num, LTR,'r',1:Num, Yzmp_SH,'b',1:Num, Yzmp,'g', 'Linewidth',lw);
% grid on
% legend('LTR','Yzmp_{SH}','Yzmp');

% Roll_BankR = u.signals.values(:,17);
% Roll_shad  = u.signals.values(:,18);
% figure
% plot(1:Num, Roll_shad,'r',1:Num, Roll_BankR,'bo', 'Linewidth',lw);
% grid on
% legend('Roll_{shad}','Roll_{BankR}');

%% Save simulation data
% SimResult_Both = u.signals.values; 
% save SimResult_Both.mat SimResult_Both
% 
% SimResult_Neither = u.signals.values; 
% save SimResult_Neither.mat SimResult_Neither
% 
% SimResult_Bank = u.signals.values; 
% save SimResult_Bank.mat SimResult_Bank
% 
% SimResult_Curvature = u.signals.values; 
% save SimResult_Curvature.mat SimResult_Curvature

