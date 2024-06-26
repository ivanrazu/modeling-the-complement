clear

% Define parameters

Kac =0.54;      
Acmax = 1e5;    
dacm = 0.05*5;     
dacn = 0.4*1.3;   
Kah = 0.6;
Ahmax = 1e5;
dahn =5;
Kc5a = 0.09;
Aifstar=3.5;
muc5a=0.1;

Kmm=0.0093;

mum=0.5184;

Knn=0.0464;

mun =1.0368;

Krn =0.05;
Krm =0.05;

mur=0.96;
sai=0.05;

Kain=0.02;
Kaim=0.04;

muai=0.15;

Kdr=0.004;
Kdn=0.0320;
mud =1.2;

dn=0.2;
dm = 0.04;
Kai=4;

drm = 2e-5;
drn=2e-5;

Kd=4;
alpha=2.5;
beta=0.5;

Kna= 0.04;
Kma = 0.04;

Kn = 3.6480;
Km = 5.4;
Kr = 45;
Kns=0.8;

Knd = 0.04;
Kmd = 0.05;

Kaid = 0.05;

% mun=0.1155;
% mum=0.0289;

gamma= 0.5;
C5astar=0.25;


Ac0=1e0*1.2;
Ah0=1;
C5a0=0;
N0=1e1*0;
M0=1e1*0;
R0=0;
Aif0=sai/muai;
D0=0;

params = [Kac, Acmax,dacm,dacn,Kah,Ahmax,dahn,Kc5a,Aifstar,muc5a,...
          Kmm,mum,Knn,mun,Krn,Krm,mur,sai,Kain,Kaim,...
          muai,Kdr,Kdn,mud,dn,dm,Kai,drm,drn,Kd,...
          alpha,beta,Kna,Kma,Kn,Km,Kr,Kns,Knd, Kmd,...
          Kaid,gamma,C5astar];

initial_conditions = [Ac0,Ah0,C5a0,N0,M0,R0,Aif0,D0];

% Time span
tspan = [0, 50];

% Call solver
sol = ode15s(@aspergillus_immune_model_v1, tspan, initial_conditions, [], params);

% Time points for which to evaluate the solution
t_eval = (0:0.1:tspan(2));

% Evaluate solution
Ys = deval(sol, t_eval);

% Define state variables
Ac = Ys(1,:);
Ah = Ys(2,:);
C5a = Ys(3,:);
N = Ys(4,:);
M = Ys(5,:);
R = Ys(6,:);
Aif = Ys(7,:);
D = Ys(8,:);

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
plot(t_eval,Ah, 'LineWidth', 2, 'DisplayName', 'Ah','LineStyle',linestyle,'Color',colors_vec{2});
% xlim([0,15])
% ylim([0,4])
xlabel('Time');
ylabel('Concentration');
legend('Location', 'northeast');
% ylim([0,500])
grid on;
%%
subplot(rows,cols,2)
hold on
plot(t_eval, C5a, 'LineWidth', 2, 'DisplayName', 'C5a','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% ylim([0,3e2])
xlim([0,tspan(end)])
grid on;

%%
subplot(rows,cols,3)
hold on
plot(t_eval, M, 'LineWidth', 2, 'DisplayName', 'M','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, N, 'LineWidth', 2, 'DisplayName', 'N','LineStyle',linestyle,'Color',colors_vec{2});
hold on
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% ylim([0,10])
xlim([0,tspan(end)])
grid on;
%%
subplot(rows,cols,4)
hold on
plot(t_eval, R, 'LineWidth', 2, 'DisplayName', 'R','LineStyle',linestyle,'Color',colors_vec{1});
hold on
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% ylim([0,40])
xlim([0,tspan(end)])
grid on;
%%
subplot(rows,cols,5)
hold on
plot(t_eval, Aif, 'LineWidth', 2, 'DisplayName', 'Aif','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% ylim([0,5])
xlim([0,tspan(end)])
grid on;
%%
subplot(rows,cols,6)
hold on
plot(t_eval, D, 'LineWidth', 2, 'DisplayName', 'D','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% ylim([0,1])
xlim([0,tspan(end)])
grid on;
%%
% return
% % Plots phagocytosis flux
return
F_aif = 1./(1+(Aif/Aifstar).^2) ;
phagoc_M=dacm * F_aif  .* Ac .* M .* R./(1+alpha*Ac);

phagoc_N = dacn * F_aif .* Ac .* N .* R./(1+alpha*Ac);

hyphae_killed_by_N =  dahn  * Ah .* N .* R .* F_aif./(1+beta*Ah);

figure
plot(t_eval,phagoc_M,'LineWidth', 2, 'DisplayName', 'Phagoc M','LineStyle',linestyle,'Color',colors_vec{1})
hold on
plot(t_eval,phagoc_N,'LineWidth', 2, 'DisplayName', 'Phagoc M','LineStyle',linestyle,'Color',colors_vec{2})
hold on
plot(t_eval,hyphae_killed_by_N,'LineWidth', 2, 'DisplayName', 'Hyphae killed by N','LineStyle',linestyle,'Color',colors_vec{3})

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

