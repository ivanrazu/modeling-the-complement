clear

% Define parameters
   
% dacm = 0.430*2;     % according to Phillippe et al 2003
dacn = 2e-1  
% Ks=1e-2*14*1e-2;
% dasm=0.4263*2; % according to Ewald 2021   % 0.1099 according to Phillippe et al 2003
% dasn=0.2354*2; % according to Ewald 2021
% Kah = 0.6931*0.1; % according to Ewald 2021
% dah=0.7895*0.1;  % according to Ewald 2021

Kc = 45*40;
Kca = 1.6e-7; % lowering Kca shortens time of spike
Kch=0.015*0.1*0.1*0.1;
muc5a = 0.24;
Kn = 17510.4;
Knn = 2.09e-08;
Kna = 0.04224;
Knd = 0.00017688;
mun = 0.0594;
Km = 12960;
Kmm = 3.72e-08;
Kma = 0.199;
Kmd = 0.00025;
mum = 0.0798;
Kd = 340;
Kdn = 1.8e-6;
Kdh = 1.4e-5;
mud = 0.228;
Kh = 40;
Khh = 1e-4;
Khd = 3e-4;
muh = 0.42;
Kcn = 1e-6*8;
Kcm = 1e-6*8;
KH = 1e-2;
Ks = 0.08;
dasn = 3e-6;
dasm = 5e-8;
Kah = 0.05;
dah = 1e-3;
Kas = 52000;
Kac = 375000;
Kahs = 0.5;
Khs = 25;
rh = 0.05;
Ahmax = 1e5;
alpha = 1e2;


Ac0=1e7;
As0=0;
Ah0=0;
C5a0=0;
N0=0;
M0=0;
D0=0;
H0=0;


M_bs=10^(4.96); %accroding to Amber's data
N_bs=10^(5.4); %%accroding to Amber's data


% M_bs=0.3*10^6 ; %accroding to Tanaka et al 2015
% N_bs=150*10^6 ; %%accroding to Tanaka et al 2015


Kn=Kn*0.1;


params = [Ks,dasm,dasn,alpha,dah,Kc,Kca,Kch,muc5a,Kn,...
         Knn,Kna,Knd,mun,Km,Kmm,Kma,Kmd,mum,Kd,...
         Kdn,Kdh,mud,Kh,Khh,Khd,muh,Kcn,Kcm,KH,...
         Kas, Kac, Kahs, Khs, rh, Ahmax];

initial_conditions = [Ac0,As0,Ah0,C5a0,N0,M0,D0,H0];

% Time span
tspan = [0, 500];

% Call solver
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
sol = ode23s(@aspergillus_immune_model_v8a, tspan, initial_conditions,[], params);

% Time points for which to evaluate the solution
t_eval = (0:0.005:tspan(2));

% Evaluate solution
Ys = deval(sol, t_eval);

% Define state variables
Ac = Ys(1,:);
As = Ys(2,:);
Ah = Ys(3,:);
C5a = Ys(4,:);
N = Ys(5,:);
M = Ys(6,:);
D = Ys(7,:);
H = Ys(8,:);
%%
% Plot the results
cols=3;
rows=2;

linestyle='-';
linewidth = 4;

colors_vec={"#000000","#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"};
%%
figure
% xSize = cols*10; X=xSize; ySize = rows*8;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
% set(gcf,'Position',[100 25 xSize*50 ySize*55]);
% 
% subplot(rows,cols,1)
% hold on
% % plot(t_eval,Ac, 'LineWidth', 2, 'DisplayName', 'Ac','LineStyle',linestyle,'Color',colors_vec{1});
% hold on;
% plot(t_eval,As, 'LineWidth', linewidth, 'DisplayName', 'As','LineStyle',linestyle,'Color',colors_vec{2});
% hold on
% plot(t_eval,Ah, 'LineWidth', linewidth, 'DisplayName', 'Ah','LineStyle',linestyle,'Color',colors_vec{3});
% xlim([0,tspan(end)])
% xlabel('Time');
% ylabel('Concentration');
% legend('Location', 'northeast');
% ylim([0,6e6])
% % grid on;
% set(gca,'linewidth',linewidth-2)
% box on;
% % yline(Kac,'--')
%%
% figure
subplot(rows,cols,1)
hold on
plot(t_eval,log10(Ac), 'LineWidth', linewidth, 'DisplayName', 'Ac','LineStyle',linestyle,'Color',colors_vec{1});
hold on;
plot(t_eval,log10(As), 'LineWidth', linewidth, 'DisplayName', 'As','LineStyle',linestyle,'Color',colors_vec{2});
hold on
plot(t_eval,log10(Ah), 'LineWidth', linewidth, 'DisplayName', 'Ah','LineStyle',linestyle,'Color',colors_vec{3});
xlim([0,tspan(end)])
xlabel('Time');
ylabel('log10(Concentration)');
legend('Location', 'Best');
ylim([-0.5,7])
% grid on;
set(gca,'linewidth',linewidth-2)
box on;

subplot(rows,cols,2)
hold on
plot(t_eval, log10(M+M_bs), 'LineWidth', linewidth, 'DisplayName', 'M','LineStyle',linestyle,'Color',colors_vec{1});
hold on
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([log10(M_bs),6])
xlim([0,tspan(end)])
ylabel('log10(Concentration)');
set(gca,'linewidth',linewidth-2)
box on;

subplot(rows,cols,3)
hold on
plot(t_eval, log10(N+N_bs), 'LineWidth', linewidth, 'DisplayName', 'N','LineStyle',linestyle,'Color',colors_vec{1});
hold on
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([log10(N_bs),6])
xlim([0,tspan(end)])
ylabel('log10(Concentration)');
set(gca,'linewidth',linewidth-2)
box on;
fontsize(16,"points")
set(gca,'linewidth',linewidth-2)

subplot(rows,cols,4)
hold on
plot(t_eval, C5a, 'LineWidth', linewidth, 'DisplayName', 'C5a','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([0,45])
xlim([0,tspan(end)])
% grid on;
set(gca,'linewidth',linewidth-2)
box on;


subplot(rows,cols,5)
hold on
plot(t_eval, D, 'LineWidth', linewidth, 'DisplayName', 'D','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([0,1.2e3])
xlim([0,tspan(end)])
% grid on;
box on;
fontsize(16,"points")
set(gca,'linewidth',linewidth-2)

subplot(rows,cols,6)
hold on
plot(t_eval, H, 'LineWidth', linewidth, 'DisplayName', 'H','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([0,1e2])
xlim([0,tspan(end)])
% grid on;
box on;
fontsize(16,"points")
set(gca,'linewidth',linewidth-2)



%%
% return
figure
xSize = 2*12; X=xSize; ySize = 1*8;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[100 25 xSize*50 ySize*55]);

subplot(1,2,1)
Ax=(0:1:Ac0);
% Kac=0.1*Kac;


f_Ax=Kas*Ax./(Kac+Ax);
semilogx(Ax,f_Ax,'LineWidth', 2)
xline(Kac,'r--','LineWidth', 2)
legend('Kas Ac /(Kac + Ac)','Ac=Kac')
legend('Location', 'Best');
xticks([0,10,100,1e3,1e4,1e5, 1e6,1e7])
xticklabels({'0','10^1','10^2','10^3','10^4','10^5','10^6','10^7'})
title('Rate of As')
xlim([0,Ac0])

subplot(1,2,2)
Ax2=(0:1:Ac0);
f_Ax2=Kahs*Ax2./(Khs+Ax2);
semilogx(Ax,f_Ax2,'LineWidth', 2, 'DisplayName','Rate of Ah')
xline(Khs,'r--','LineWidth', 2)
legend('Kahs As /(Khs + As)','As=Khs')
legend('Location', 'Best');
fontsize(16,"points")
xticks([0,10,100,1e3,1e4,1e5, 1e6,1e7])
xticklabels({'0','10^1','10^2','10^3','10^4','10^5','10^6','10^7'})
title('Rate of Ah')

%%



% 
% subplot(rows,cols,6)
% hold on
% plot(t_eval,log10(Ac), 'LineWidth', 2, 'DisplayName', 'Ac','LineStyle',linestyle,'Color',colors_vec{1});
% hold on;
% plot(t_eval,log10(As), 'LineWidth', 2, 'DisplayName', 'As','LineStyle',linestyle,'Color',colors_vec{2});
% hold on
% plot(t_eval,log10(Ah), 'LineWidth', 2, 'DisplayName', 'Ah','LineStyle',linestyle,'Color',colors_vec{3});
% 
% xlim([0,tspan(end)])
% xlabel('Time');
% ylabel('Concentration');
% legend('Location', 'northeast');
% % ylim([0,1e6])
% grid on;



%%
% return
% % Plots phagocytosis of Ac flux
% F_aif = 1./(1+(Aif/Aifstar).^2) ;
% phagoc_Ac_by_M = dacm * F_aif  .* Ac .* M .* R./(1+alpha*Ac);
% 
% phagoc_Ac_by_N = dacn * F_aif .* Ac .* N .* R./(1+alpha*Ac);
% 
% 
% phagoc_As_by_M = dasm * F_aif  .* As .* M .* R./(1+beta*As);
% 
% phagoc_As_by_N = dasn * F_aif .* As .* N .* R./(1+beta*As);
% 
% 
% hyphae_killed_by_N =  dah  * Ah .* N .* R .* F_aif./(1+gamma*Ah);
% 
% figure
% plot(t_eval,phagoc_Ac_by_M,'LineWidth', 2, 'DisplayName', 'Phagoc Ac by M','LineStyle',linestyle,'Color',colors_vec{1})
% hold on
% plot(t_eval,phagoc_Ac_by_N,'LineWidth', 2, 'DisplayName', 'Phagoc Ac by N','LineStyle',linestyle,'Color',colors_vec{2})
% hold on
% plot(t_eval,phagoc_As_by_N,'LineWidth', 2, 'DisplayName', 'Phagoc As by N','LineStyle',linestyle,'Color',colors_vec{3})
% hold on
% plot(t_eval,phagoc_As_by_M,'LineWidth', 2, 'DisplayName', 'Phagoc As by M','LineStyle',linestyle,'Color',colors_vec{4})
% hold on
% plot(t_eval,hyphae_killed_by_N,'LineWidth', 2, 'DisplayName', 'Hyphae killed by N','LineStyle',linestyle,'Color',colors_vec{5})
% 
% xlabel('Time');
% legend('Location', 'Best');
% grid on;




%%
% return

% Ndat = xlsread('Neutrophils_WT.xlsx','Sheet1', 'A14:E18');
% meanNdat = mean(Ndat);
% sdNdat = std(Ndat);
% 
% 
% Mdat = xlsread('Monocytes_WT.xlsx','Sheet1', 'A12:E16');
% meanMdat = mean(Mdat);
% sdMdat = std(Mdat);
% 
% CXCL2_dat = xlsread('CXCL2_WT.xlsx','Sheet1', 'A12:E16');
% meanCXCL2dat = mean(CXCL2_dat);
% sdCXCL2dat = std(CXCL2_dat);
% 
% tscat = [0,12,24,48,72];
% 
% 
% figure
% 
% subplot(1,3,1)
% errorbar(tscat,meanNdat,sdNdat,'ks','MarkerSize',10,'MarkerFaceColor','black','LineStyle','none','LineWidth',2);
% hold on;
% xlabel('Time (d)');
% ylabel('Log_{10} N');
% xticks([0,12,24,48,72])
% xticklabels({'0','12','24','48','72'})
% grid on
% % xlim([-0.2,12])
% % ylim([3.5,5.5])
% xtickangle(0)
% hold on
% 
% 
% subplot(1,3,2)
% errorbar(tscat,meanCXCL2dat,sdCXCL2dat,'ks','MarkerSize',10,'MarkerFaceColor','black','LineStyle','none','LineWidth',2);
% hold on;
% xlabel('Time (d)');
% ylabel('Log_{10} CXXL2');
% xticks([0,12,24,48,72])
% xticklabels({'0','12','24','48','72'})
% grid on
% % xlim([-0.2,12])
% % ylim([3.5,5.5])
% xtickangle(0)
% hold on
% 
% subplot(1,3,3)
% errorbar(tscat,meanMdat,sdMdat,'ks','MarkerSize',10,'MarkerFaceColor','black','LineStyle','none','LineWidth',2);
% hold on;
% xlabel('Time (d)');
% ylabel('Log_{10} M');
% xticks([0,12,24,48,72])
% xticklabels({'0','12','24','48','72'})
% grid on
% % xlim([-0.2,12])
% % ylim([3.5,5.5])
% xtickangle(0)
% hold on
%%
% x = (0:0.1:7);
% xstar = 3.5;
% 
% fx= 1./(1+(x/xstar).^4);
% 
% 
% figure
% semilogx(x,fx)
% hold on
% xline(xstar,'r--')
%%

% Nx=(0:0.1:30);
% Mx=(0:0.1:30);
% 
% % f_ROS_N=Kr *(Krn*N)./(1+(Krn*N));
% f_ROS_N=Kr *(Krn*Nx)./(1+(Krn*Nx));
% % f_ROS_M = Kr *(Krm*M)./(1+(Krm*M));
% f_ROS_M = Kr *(Krm*Mx)./(1+(Krm*Mx));
% 
% figure
% plot(Nx,f_ROS_N,'LineWidth', 2, 'DisplayName', 'ROS by N','LineStyle',linestyle,'Color',colors_vec{2})
% hold on
% plot(Mx,f_ROS_M,'LineWidth', 2, 'DisplayName', 'ROS by M','LineStyle','--','Color',colors_vec{1})
% xlabel('Time');
% legend('Location', 'Best');
% grid on;

