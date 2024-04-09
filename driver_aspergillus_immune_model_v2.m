clear

% Define parameters
   
dacm = 0.25*0.5;     
dacn = 0.8*0.5;  
alpha=2.5*0.5*0.5*0;
Ks=0.03*10*5*0.5;
dasm=0.038;%0.01*0.5;
dasn=0.05*0.5;
beta=0.5*0;
Kah = 0.1*10*6*0.2*0.5;
dah=0.8;
gamma= 0.5*0.5*0;

Kc = 0.09*5;
Kca = 0.04*4;
Kch=0.003*5;
muc5a=0.1;

Kn = 3.6480*1.2;%*0.1;
Knn=0.4640*4;
Kna= 0.04*2*2*2;
Knd = 0.04*2*2;
mun =1.03688*0.8;


dnc =0.5;

dns=0.1*5;
dnh=0.5;
Km = 5.4*0.6;
Kmm=0.0930*10*4;
Kma = 0.04;  
Kmd = 0.05;
mum=0.5184*1.5;
dmc = 0.02;
dms=0.01;
Kr = 45*0.5;

Krn =0.05;
Krm =0.05;
mur=0.96;
drm = 2e-5;
drn=2e-5;
sai=0.05;
Kai=4;
Kain=0.02;
Kaim=0.04;
Kaid = 0.05;

muai=0.15;
Kd=4;
Kdr=0.004;
Kdn=0.0320;
Kdh = 0.02;
mud =1.2;
Kh = 0.1*4;
Khh=0.02*4;
Khd = 0.03*4;
muh=0.01*15;

Aifstar=3.5;
C5astar=0.25;


Ac0=7;
As0=0;
Ah0=0;
C5a0=0;
N0=1e1*0;
M0=1e1*0;
R0=1;
Aif0=sai/muai;
D0=0;
H0=0;

params = [dacm,dacn,alpha,Ks,dasm,dasn,beta,Kah,dah,gamma,...
          Kc,Kca,Kch,muc5a,Kn,Knn,Kna,Knd,mun,dnc...
          dns,dnh,Km,Kmm,Kma,Kmd,mum,dmc,dms,Kr,...
          Krn,Krm,mur,drm,drn,sai,Kai,Kain,Kaim,Kaid,...
          muai,Kd,Kdr,Kdn,Kdh,mud,Kh,Khh,Khd,muh,...
          Aifstar,C5astar];

initial_conditions = [Ac0,As0,Ah0,C5a0,N0,M0,R0,Aif0,D0,H0];

% Time span
tspan = [0, 20];

% Call solver
sol = ode15s(@aspergillus_immune_model_v2, tspan, initial_conditions, [], params);

% Time points for which to evaluate the solution
t_eval = (0:0.1:tspan(2));

% Evaluate solution
Ys = deval(sol, t_eval);

% Define state variables
Ac = Ys(1,:);
As = Ys(2,:);
Ah = Ys(3,:);
C5a = Ys(4,:);
N = Ys(5,:);
M = Ys(6,:);
R = Ys(7,:);
Aif = Ys(8,:);
D = Ys(9,:);
H = Ys(10,:);
%%
% Plot the results
cols=3;
rows=2;

linestyle='-';

colors_vec={"#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"};

figure
subplot(rows,cols,1)
hold on
plot(t_eval,Ac, 'LineWidth', 2, 'DisplayName', 'Ac','LineStyle',linestyle,'Color',colors_vec{1});
hold on;
plot(t_eval,As, 'LineWidth', 2, 'DisplayName', 'As','LineStyle',linestyle,'Color',colors_vec{2});
hold on
plot(t_eval,Ah, 'LineWidth', 2, 'DisplayName', 'Ah','LineStyle',linestyle,'Color',colors_vec{3});

% xlim([0,15])
ylim([0,7.5])
xlabel('Time');
ylabel('Concentration');
legend('Location', 'northeast');
% ylim([0,500])
grid on;

subplot(rows,cols,2)
hold on
plot(t_eval, C5a, 'LineWidth', 2, 'DisplayName', 'C5a','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([0,0.2])
xlim([0,tspan(end)])
grid on;


subplot(rows,cols,3)
hold on
plot(t_eval, M, 'LineWidth', 2, 'DisplayName', 'M','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, N, 'LineWidth', 2, 'DisplayName', 'N','LineStyle',linestyle,'Color',colors_vec{2});
hold on
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([0,2.5])
xlim([0,tspan(end)])
grid on;

subplot(rows,cols,4)
hold on
plot(t_eval, R, 'LineWidth', 2, 'DisplayName', 'R','LineStyle',linestyle,'Color',colors_vec{1});
hold on
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([0,2.5])
xlim([0,tspan(end)])
grid on;

subplot(rows,cols,5)
hold on
plot(t_eval, Aif, 'LineWidth', 2, 'DisplayName', 'Aif','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([0,2])
xlim([0,tspan(end)])
grid on;

subplot(rows,cols,6)
hold on
plot(t_eval, D, 'LineWidth', 2, 'DisplayName', 'D','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, H, 'LineWidth', 2, 'DisplayName', 'H','LineStyle',linestyle,'Color',colors_vec{2});

xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([0,0.15])
xlim([0,tspan(end)])
grid on;

return

%%
% return
% % Plots phagocytosis of Ac flux
F_aif = 1./(1+(Aif/Aifstar).^2) ;
phagoc_Ac_by_M = dacm * F_aif  .* Ac .* M .* R./(1+alpha*Ac);

phagoc_Ac_by_N = dacn * F_aif .* Ac .* N .* R./(1+alpha*Ac);


phagoc_As_by_M = dasm * F_aif  .* As .* M .* R./(1+beta*As);

phagoc_As_by_N = dasn * F_aif .* As .* N .* R./(1+beta*As);


hyphae_killed_by_N =  dah  * Ah .* N .* R .* F_aif./(1+gamma*Ah);

figure
plot(t_eval,phagoc_Ac_by_M,'LineWidth', 2, 'DisplayName', 'Phagoc Ac by M','LineStyle',linestyle,'Color',colors_vec{1})
hold on
plot(t_eval,phagoc_Ac_by_N,'LineWidth', 2, 'DisplayName', 'Phagoc Ac by N','LineStyle',linestyle,'Color',colors_vec{2})
hold on
plot(t_eval,phagoc_As_by_N,'LineWidth', 2, 'DisplayName', 'Phagoc As by N','LineStyle',linestyle,'Color',colors_vec{3})
hold on
plot(t_eval,phagoc_As_by_M,'LineWidth', 2, 'DisplayName', 'Phagoc As by M','LineStyle',linestyle,'Color',colors_vec{4})
hold on
plot(t_eval,hyphae_killed_by_N,'LineWidth', 2, 'DisplayName', 'Hyphae killed by N','LineStyle',linestyle,'Color',colors_vec{5})

xlabel('Time');
legend('Location', 'Best');
grid on;




%%
% return

Ndat = xlsread('Neutrophils_WT.xlsx','Sheet1', 'A14:E18');
meanNdat = mean(Ndat);
sdNdat = std(Ndat);


Mdat = xlsread('Monocytes_WT.xlsx','Sheet1', 'A12:E16');
meanMdat = mean(Mdat);
sdMdat = std(Mdat);

CXCL2_dat = xlsread('CXCL2_WT.xlsx','Sheet1', 'A12:E16');
meanCXCL2dat = mean(CXCL2_dat);
sdCXCL2dat = std(CXCL2_dat);

tscat = [0,12,24,48,72];


figure

subplot(1,3,1)
errorbar(tscat,meanNdat,sdNdat,'ks','MarkerSize',10,'MarkerFaceColor','black','LineStyle','none','LineWidth',2);
hold on;
xlabel('Time (d)');
ylabel('Log_{10} N');
xticks([0,12,24,48,72])
xticklabels({'0','12','24','48','72'})
grid on
% xlim([-0.2,12])
% ylim([3.5,5.5])
xtickangle(0)
hold on


subplot(1,3,2)
errorbar(tscat,meanCXCL2dat,sdCXCL2dat,'ks','MarkerSize',10,'MarkerFaceColor','black','LineStyle','none','LineWidth',2);
hold on;
xlabel('Time (d)');
ylabel('Log_{10} CXXL2');
xticks([0,12,24,48,72])
xticklabels({'0','12','24','48','72'})
grid on
% xlim([-0.2,12])
% ylim([3.5,5.5])
xtickangle(0)
hold on

subplot(1,3,3)
errorbar(tscat,meanMdat,sdMdat,'ks','MarkerSize',10,'MarkerFaceColor','black','LineStyle','none','LineWidth',2);
hold on;
xlabel('Time (d)');
ylabel('Log_{10} M');
xticks([0,12,24,48,72])
xticklabels({'0','12','24','48','72'})
grid on
% xlim([-0.2,12])
% ylim([3.5,5.5])
xtickangle(0)
hold on
%%
x = (0:0.1:7);
xstar = 3.5;

fx= 1./(1+(x/xstar).^4);


figure
semilogx(x,fx)
hold on
xline(xstar,'r--')
%%

Nx=(0:0.1:30);
Mx=(0:0.1:30);

% f_ROS_N=Kr *(Krn*N)./(1+(Krn*N));
f_ROS_N=Kr *(Krn*Nx)./(1+(Krn*Nx));
% f_ROS_M = Kr *(Krm*M)./(1+(Krm*M));
f_ROS_M = Kr *(Krm*Mx)./(1+(Krm*Mx));

figure
plot(Nx,f_ROS_N,'LineWidth', 2, 'DisplayName', 'ROS by N','LineStyle',linestyle,'Color',colors_vec{2})
hold on
plot(Mx,f_ROS_M,'LineWidth', 2, 'DisplayName', 'ROS by M','LineStyle','--','Color',colors_vec{1})
xlabel('Time');
legend('Location', 'Best');
grid on;



