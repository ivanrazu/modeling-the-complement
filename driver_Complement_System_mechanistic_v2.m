clear

% Initial conditions
C1_0 = 45 * 10^(-2);
C2_0 = 21.36 * 10^(-2);
C3_0 = 888.9 * 10^(-2);
C4_0 = 310.67 * 10^(-2);
C5_0 = 44.44 * 10^(-2);

C3convCPLP_0 = 0;
C3convAP_0 = 0;

C5conv_0 = 0;
C3a_0 = 0;
C5a_0 = 0;
C3b_0 = 0;
C3H2O_0 = 0;

Ac_0=1e2*1 ;
Ah_0 = 1e1*1;

Am_0 = 0;
N_0 = 0;

C1inh_0 = 1.7; % (uM according to Korotaevskiy)
FB_0 = 2.1; % Factor B (uM according to Korotaevskiy)
FH_0 = 3.2; % Factor H (uM according to Korotaevskiy)
FI_0 = 0.4; % Factor I (uM according to Korotaevskiy)
P_0 = 0.42; % Properdin (uM according to Korotaevskiy)


% Define parameters

Kc1 = 0.00021;
Kc2 = 0.0002;
Kc3 = 0.0005;
Kc4 = 0.0002;
Kc5 = 0.0002;

muc1 = 1e-1*1;
muc2 = 1e-1*1;
muc3 = 1e-1*1;
muc4 = 1e-1*1;
muc5 = 1e-1*1;

sc1 = C1_0*muc1;
sc2 = C2_0*muc2;
sc3 = C3_0*muc3;
sc4 = C4_0*muc4;
sc5 = C5_0*muc5;

C1star=1e0;
C2star=1e0;
C3star=1e0;
C4star=1e0;
C5star=1e0;

Kc3convcpcl = 1e-6;
Kcp = 0.01;
Klp= 0.2;

C1inhstar = 1;
muc3convcplp = 1e-2;

Kc3convap = 1e-1*4;
C3H2Ostar= 1e0;
muc3convap = 3e-2;

Kc5conv = 1e-1;
Kc5convhs = 1e1;
muc5conv = 1e-2;

Kc3acplp = 0.05;
Kc3aap = 0.02;
muc3a = 6e-2;

Kc5a = 5e-2;
muc5a = 1e-2;

Kc3bcplp = 0.05;
Kc3bh2ofb = 0.03;
Kc3bap = 0.03;

FHstar = 1;
FIstar = 1;
muc3b = 4e-2;

Kc3h2o = 0.0001; % 1/min  K1 in Bakshi et al 2020, reference Pangburn et al 1981
muc3h2o = 5e-2;

KAc = 0.1;
Acmax = 1e4;

dacm = 5e-2;
dacn = 1e-3;

KAh = 0.25;
Ahmax = 1e4;

dahm = 6e-3;
dahn = 1e-2;

Kam = 0.1;
Kamc3a = 0.5;
Kamc5a = 0.5;
Amhs = 1e2;
muam = 1e-2*5;

dam = 1e-3;
Kn = 0.2;
Knc3a = 0.5;
Knc5a = 0.5;
Nhs = 1e2;
mun = 1e-2;
dan = 1e-1;

sc1inh = 0.2*0 ;
muc1inh = 1e-1*0 ;

sfb = 0.000798; % muM/min, ks2 in Bakshi et al, reference Alper and Rosen 1984
mufb = 0.000333; % 1/min d2 in Bakshi et al, reference Alper and Rosen 1984
sfh = 0.00067;  % muM/min, ks3 in Bakshi et al, reference Charlesworth et al 1979, Dopler et al 2019
mufh = 0.00022; % 1/min d3 in Bakshi et al, reference Charlesworth et al 1979, Dopler et al 2019
sfi = 0.2 *0;
mufi = 1e-2 *0;

Kp=0.00007; % muM/min ks4 in Bakshi et al, reference Ziegler et al 1975
mup=0.000016; % 1/min d4 in Bakshi et al, reference Ziegler et al 1975
Pstar=1e0;



params = [Kc1, Kc2, Kc3, Kc4, Kc5, sc1, sc2, sc3, sc4, sc5, muc1, muc2, muc3, muc4, muc5,...
    C1star, C2star, C3star, C4star, C5star, Kc3convcpcl, Kcp, Klp, C1inhstar,...
    muc3convcplp, Kc3convap, C3H2Ostar, muc3convap, Kc5conv, Kc5convhs, muc5conv,...
    Kc3acplp, Kc3aap, muc3a, Kc5a, muc5a, Kc3bcplp, Kc3bh2ofb, Kc3bap, FHstar,...
    FIstar, muc3b, Kc3h2o, muc3h2o, KAc, Acmax, dacm, dacn, KAh, Ahmax, dahm, dahn,...
    Kam, Kamc3a, Kamc5a, Amhs, muam, dam, Kn, Knc3a, Knc5a, Nhs, mun, dan,...
    sc1inh, muc1inh, sfb, mufb, sfh, mufh, sfi, mufi,Kp,mup,Pstar];

initial_conditions = [C1_0, C2_0, C3_0, C4_0, C5_0, ...
  C3convCPLP_0, C3convAP_0, C5conv_0, C3a_0,C5a_0,C3b_0, C3H2O_0, Ac_0, ...
   Ah_0,Am_0, N_0, C1inh_0, FB_0, FH_0, FI_0,P_0];

% Time span
tspan = [0, 250];

% Call solver
sol = ode23s(@Complement_System_mechanistic_v2, tspan, initial_conditions, [], params);

% Time points for which to evaluate the solution
t_eval = (0:0.001:tspan(2));

% Evaluate solution
Ys = deval(sol, t_eval);

%% Define state variables frrom solutions

C1 = Ys(1, :);
C2 = Ys(2, :);
C3 = Ys(3, :);
C4 = Ys(4, :);
C5 = Ys(5, :);
C3convCPLP = Ys(6, :);
C3convAP =Ys(7,:);
C5conv = Ys(8, :);
C3a = Ys(9, :);
C5a = Ys(10, :);
C3b = Ys(11, :);
C3H2O = Ys(12, :);
Ac = Ys(13, :);
Ah = Ys(14,:);
Am = Ys(15,:);
N = Ys(16,:);
C1inh =Ys(17,:) ;
FB = Ys(18,:); % Factor B
FH = Ys(19,:); % Factor H
FI = Ys(20,:); % Factor I
P = Ys(21,:); % Properdin

%%
% Plot the results
cols=3;
rows=2;

linestyle='-';

colors_vec={"#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"};

figure
subplot(rows,cols,1)
hold on
plot(t_eval,C1, 'LineWidth', 2, 'DisplayName', 'C1','LineStyle',linestyle,'Color',colors_vec{1});
hold on;
plot(t_eval, C2, 'LineWidth', 2, 'DisplayName', 'C2','LineStyle',linestyle,'Color',colors_vec{3});
hold on
plot(t_eval, C3, 'LineWidth', 2, 'DisplayName', 'C3','LineStyle',linestyle,'Color',colors_vec{4});
hold on
plot(t_eval, C4, 'LineWidth', 2, 'DisplayName', 'C4','LineStyle',linestyle,'Color',colors_vec{5});
hold on
plot(t_eval, C5, 'LineWidth', 2, 'DisplayName', 'C5','LineStyle',linestyle,'Color',colors_vec{6});
hold on
plot(t_eval, C3H2O, 'LineWidth', 2, 'DisplayName', 'C3H2O','LineStyle',linestyle,'Color',colors_vec{7});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'northeast');
% title('Classic Pathway Dynamics');
grid on;

%%
subplot(rows,cols,2)
hold on
plot(t_eval, C3convCPLP, 'LineWidth', 2, 'DisplayName', 'C3convCPLP','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, C3convAP, 'LineWidth', 2, 'DisplayName', 'C3convAP','LineStyle',linestyle,'Color',colors_vec{2});
hold on
plot(t_eval, C5conv, 'LineWidth', 2, 'DisplayName', 'C5conv','LineStyle',linestyle,'Color',colors_vec{3});

xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% title('Classic Pathway Dynamics');
grid on;
% ylim([0,10e-12])

%%
subplot(rows,cols,3)
hold on
plot(t_eval, Ac, 'LineWidth', 2, 'DisplayName', 'Asperg. conidia','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, Ah, 'LineWidth', 2, 'DisplayName', 'Asperg. hyphae','LineStyle',linestyle,'Color',colors_vec{2});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% title('Classic Pathway Dynamics');
grid on;

%%
subplot(rows,cols,4)
plot(t_eval, C3a, 'LineWidth', 2, 'DisplayName', 'C3a','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval,C5a, 'LineWidth', 2, 'DisplayName', 'C5a','LineStyle',linestyle,'Color',colors_vec{2});
hold on
plot(t_eval,C3b, 'LineWidth', 2, 'DisplayName', 'C3b','LineStyle',linestyle,'Color',colors_vec{3});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% title('Classic Pathway Dynamics');
grid on;

%%

subplot(rows,cols,5)
hold on
plot(t_eval, C1inh, 'LineWidth', 2, 'DisplayName', 'C1inh','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, FB, 'LineWidth', 2, 'DisplayName', 'FB','LineStyle',linestyle,'Color',colors_vec{2});
hold on
plot(t_eval, FH, 'LineWidth', 2, 'DisplayName', 'FH','LineStyle',linestyle,'Color',colors_vec{3});
hold on
plot(t_eval, FI, 'LineWidth', 2, 'DisplayName', 'FI','LineStyle',linestyle,'Color',colors_vec{4});
hold on
plot(t_eval, P, 'LineWidth', 2, 'DisplayName', 'P','LineStyle',linestyle,'Color',colors_vec{5});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% title('Classic Pathway Dynamics');
grid on;
hold off;
%%

subplot(rows,cols,6)
hold on
plot(t_eval, Am, 'LineWidth', 2, 'DisplayName', 'Am','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, N, 'LineWidth', 2, 'DisplayName', 'N','LineStyle',linestyle,'Color',colors_vec{2});
hold on
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
% title('Classic Pathway Dynamics');
grid on;


return

%% 
figure
subplot(rows,cols,1)
hold on
% plot(t_eval, C4, 'LineWidth', 2, 'DisplayName', 'C4','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, C1, 'LineWidth', 2, 'DisplayName', 'C1','LineStyle',linestyle,'Color',colors_vec{2});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
%%
subplot(rows,cols,2)
hold on
plot(t_eval, AgAb, 'LineWidth', 2, 'DisplayName', 'AgAb','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;

subplot(rows,cols,3)
hold on
plot(t_eval, C4bC2aC3bC5, 'LineWidth', 2, 'DisplayName', 'C4bC2aC3bC5','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;

subplot(rows,cols,4)
hold on
plot(t_eval, C1bar, 'LineWidth', 2, 'DisplayName', 'C1bar','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Concentration');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;




