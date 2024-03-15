clear

% Define parameters

Kac = 0.0250*4;      
Acmax = 1e5;    
dacm = 0.012*5*5;     
dacn = 0.06;   

Kah = 0.03*1.8;
Ahmax = 1e5;
dahn =0.06;
Kc5a = 0.1;
Aifstar=1e1;
muc5a=0.1;

Km=0.0320;

mum=0.0224*2*2*2;

Kn=0.0464;

mun =0.0324*2*2*2;

Krn =0.05;
Krm =0.05;

mur=0.06*13;
sai=0.3;
Kain=0.02;
Kaim=0.04;
muai=0.03;
Kdr=0.004;
Kdn=0.0320;
mud =1.2;

dn=0.02;
dm = 0.02;
Kai=4;

drm = 0.02*0.05;
drn=0.02*0.05;

Kd=4;
alpha=2.5;
beta=1.5;

Kna= 0.05;
Kma = 0.05;

Ac0=1e2;
Ah0=0;
C5a0=0;
N0=1e1*0;
M0=1e1*0;
R0=0;
Aif0=sai/muai;
D0=0;

params = [Kac, Acmax,dacm,dacn,Kah,Ahmax,dahn,Kc5a,Aifstar,muc5a,...
    Km,mum,Kn,mun,Krn,Krm,mur,sai,Kain,Kaim,muai,Kdr,Kdn,mud,dn,...
    dm,Kai,drm,drn,Kd,alpha,beta,Kna,Kma];

initial_conditions = [Ac0,Ah0,C5a0,N0,M0,R0,Aif0,D0];

% Time span
tspan = [0, 100];

% Call solver
sol = ode15s(@aspergillus_immune_model_v1, tspan, initial_conditions, [], params);

% Time points for which to evaluate the solution
t_eval = (0:0.0001:tspan(2));

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
xlim([0,20])
xlabel('Time');
ylabel('Concentration');
legend('Location', 'northeast');
grid on;
%%
subplot(rows,cols,2)
hold on
plot(t_eval, C5a, 'LineWidth', 2, 'DisplayName', 'C5a','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
grid on;

%%
subplot(rows,cols,3)
hold on
plot(t_eval, M, 'LineWidth', 2, 'DisplayName', 'M','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, N, 'LineWidth', 2, 'DisplayName', 'N','LineStyle',linestyle,'Color',colors_vec{2});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([0,100])
grid on;
%%
subplot(rows,cols,4)
hold on
plot(t_eval, R, 'LineWidth', 2, 'DisplayName', 'R','LineStyle',linestyle,'Color',colors_vec{1});
hold on
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
ylim([0,40])
grid on;
%%
subplot(rows,cols,5)
hold on
plot(t_eval, Aif, 'LineWidth', 2, 'DisplayName', 'Aif','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
grid on;
%%
subplot(rows,cols,6)
hold on
plot(t_eval, D, 'LineWidth', 2, 'DisplayName', 'D','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
grid on;

