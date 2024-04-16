clear

% Define parameters
   
dacm = 0.0430*0.001*0;     % according to Phillippe et al 2003
dacn = 1e-5*0;  
Ks=10^6;
dasm=0.4263*0.1; % according to Ewald 2021   % 0.1099 according to Phillippe et al 2003
dasn=0.2354*0.1; % according to Ewald 2021
Kah = 0.6931; % according to Ewald 2021
dah=0.7895;  % according to Ewald 2021
Kc = 0.45*1000*10;
Kca = 0.16*0.00001;
Kch=0.015;
muc5a=0.1*8*2;
Kn = 43.776*2;
Knn=0.0019;
Kna= 9.6e-5;
Knd = 1.608e-4;
mun =0.0594*10*2; % according to Ewald 2021 
dnc =0.5475*2; % according to Ewald 2021 
dns=0.5475*2;  % according to Ewald 2021 
dnh=0.5475*2;  % according to Ewald 2021 
Km = 3.24;
Kmm=3.72e-04;
Kma = 1e-5;  
Kmd = 2.5e-5;
mum=0.0798; % according to Ewald 2021 
sai=2.3232; % so that sai/muai = 10^(1.19) which is IL10-0 in Ambers' data
Kai=4;
Kain=0.02;
Kaim=0.04;
Kaid = 0.05;
muai=0.15*1.7;
Kd=1.7;
Kdn=0.0018;
Kdh = 1.4e-7;
mud =2.28;
Kh = 0.2;
Khh=1e-4;
Khd = 3e-4;
muh=0.9;
Aifstar=80;

% Kn=Kn*0.1;

Ac0=1e6;
As0=0;
Ah0=0;
C5a0=0;
N0=0;
M0=0;
Aif0=sai/muai;
D0=0;
H0=0;


% M0=10^(4.96)*0; %accroding to Amber's data
% N0=10^(5.4)*0; %%accroding to Amber's data



params = [dacm,dacn,Ks,dasm,dasn,Kah,dah,Kc,Kca,Kch,...
         muc5a,Kn,Knn,Kna,Knd,mun,Km,Kmm,Kma,Kmd,...
         mum,sai,Kai,Kain,Kaim,Kaid,muai,Kd,Kdn,Kdh,...
         mud,Kh,Khh,Khd,muh,Aifstar];

initial_conditions = [Ac0,As0,Ah0,C5a0,N0,M0,Aif0,D0,H0];

% Time span
tspan = [0,  20];

% Call solver
sol = ode15s(@aspergillus_immune_model_v6, tspan, initial_conditions, [], params);

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
Aif = Ys(7,:);
D = Ys(8,:);
H = Ys(9,:);
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

xlim([0,tspan(end)])
xlabel('Time');
ylabel('Concentration');
legend('Location', 'northeast');
ylim([0,1e6])
grid on;

subplot(rows,cols,2)
hold on
plot(t_eval, C5a, 'LineWidth', 2, 'DisplayName', 'C5a','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% ylim([0,1.5])
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
ylim([0,150])
xlim([0,tspan(end)])
grid on;

% subplot(rows,cols,4)
% hold on
% plot(t_eval, R, 'LineWidth', 2, 'DisplayName', 'R','LineStyle',linestyle,'Color',colors_vec{1});
% hold on
% xlabel('Time');
% ylabel('Concentration');
% legend('Location', 'Best');
% ylim([0,2.5])
% xlim([0,tspan(end)])
% grid on;

subplot(rows,cols,4)
hold on
plot(t_eval, Aif, 'LineWidth', 2, 'DisplayName', 'Aif','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% ylim([0,25])
xlim([0,tspan(end)])
grid on;

subplot(rows,cols,5)
hold on
plot(t_eval, D, 'LineWidth', 2, 'DisplayName', 'D','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, H, 'LineWidth', 2, 'DisplayName', 'H','LineStyle',linestyle,'Color',colors_vec{2});

xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% ylim([0,0.5])
xlim([0,tspan(end)])
grid on;


subplot(rows,cols,6)
hold on
plot(t_eval,log10(Ac), 'LineWidth', 2, 'DisplayName', 'Ac','LineStyle',linestyle,'Color',colors_vec{1});
hold on;
plot(t_eval,log10(As), 'LineWidth', 2, 'DisplayName', 'As','LineStyle',linestyle,'Color',colors_vec{2});
hold on
plot(t_eval,log10(Ah), 'LineWidth', 2, 'DisplayName', 'Ah','LineStyle',linestyle,'Color',colors_vec{3});

xlim([0,tspan(end)])
xlabel('Time');
ylabel('Concentration');
legend('Location', 'northeast');
% ylim([0,1e6])
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

%%
figure
Ax=(0:1:10000);
alphax=2000;
f_Ax=Ax./(alphax+Ax);
semilogx(Ax,f_Ax,'LineWidth', 2)
xline(alphax,'r--','LineWidth', 2)