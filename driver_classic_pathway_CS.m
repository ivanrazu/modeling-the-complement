clear

% Define parameters
% K0 = 0.5 * 10^(1); % arbitrarily chosen to obtain biologically relevant trajectory
K0=10.7722;

% Parameters taken from Hirayama et al 1996  
K1 = 0.5 * 10^8;
K_1 = 0.12 * 10^(-5);
K2 = 0.96 * 10^6;
K3 = 1.27 * 10^6;
K_3 = 0.84 * 10^(-6);
K4 = 2.04 * 10^6;
K5 = 1.6 * 10^6;
K_5 = 1.15 * 10^(-6);
K6 = 3.1 * 10^6;
K7 = 2.1 * 10^8;
K_7 = 0.43 * 10^(-6);
K8 = 4.8 * 10^8;

K9 = 8 * 10^(7); % arbitrarily chosen to obtain biologically relevant trajectory
K_9 = 0.3 * 10^(-6); % arbitrarily chosen to obtain biologically relevant trajectory
 

KAgAb = 0.5 * 10^(1);
KC1inh=0.06*1e3;
C1inh=1.9231*10^(-6); % Taken from Hirayama 

FH=2.25*10^(-6)*0 ;% Taken from Hirayama
C4bp=7.01754*10^(-8)*0;% Taken from Hirayama

DAF=0;
CR1=0;
CR2=0;
MCP=0;

% KAgAb=0;

params = [K0, K1, K_1, K2, K3, K_3, K4, K5, K_5, K6, K7, K_7, K8, K9, K_9,FH,C4bp,DAF,CR1,CR2,MCP,KAgAb,KC1inh,C1inh];

% Initial conditions
C1_0 = 45 * 10^(-5);
C1bar_0 = 0;
C4_0 = 310.67 * 10^(-5);
C1barC4_0 = 0;
C4a_0 = 0;
C4b_0 = 0;
C2_0 = 21.36 * 10^(-5);
C4b2_0 = 0;
C4bC2a_0 = 0;
C2b_0 = 0;
C3_0 = 888.9 * 10^(-5);
C4bC2aC3_0 = 0;
C4bC2ac3b_0 = 0;
C3a_0 = 0;
C5_0 = 44.44 * 10^(-5);
C4bC2aC3bC5_0 = 0;
C4bC2aC3bC5b_0 = 0;
C5a_0 = 0;
C6789_0 = 2*10^(-4);  % arbitrarily chosen to obtain biologically relevant trajectory
MAC_0 = 0;

AgAb_0=1e2;

initial_conditions = [C1_0, C1bar_0, C4_0, C1barC4_0, C4a_0, C4b_0, C2_0, ...
    C4b2_0, C4bC2a_0, C2b_0, C3_0, C4bC2aC3_0, C4bC2ac3b_0, ...
    C3a_0, C5_0, C4bC2aC3bC5_0, C4bC2aC3bC5b_0, C5a_0, C6789_0, MAC_0,AgAb_0]/1000;

% Time span
tspan = [0, 2];

% Call solver
sol = ode23s(@classic_pathway_CS, tspan, initial_conditions, [], params);

% Time points for which to evaluate the solution
t_eval = linspace(tspan(1), tspan(2), 1000);

% Evaluate solution
Ys = deval(sol, t_eval);

%% Define state variables frrom solutions

C1 = Ys(1, :);
C1bar = Ys(2, :);
C4 = Ys(3, :);
C1barC4 = Ys(4, :);
C4a = Ys(5, :);
C4b = Ys(6, :);
C2 = Ys(7, :);
C4b2 = Ys(8, :);
C4bC2a = Ys(9, :);
C2b = Ys(10, :);
C3 = Ys(11, :);
C4bC2aC3 = Ys(12, :);
C4bC2aC3b = Ys(13, :);
C3a = Ys(14, :);
C5 = Ys(15, :);
C4bC2aC3bC5 = Ys(16, :);
C4bC2aC3bC5b = Ys(17, :);
C5a = Ys(18, :);
C6C7C8C9 = Ys(19, :);
MAC = Ys(20, :);
AgAb = Ys(21, :);


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
plot(t_eval,C1bar, 'LineWidth', 2, 'DisplayName', 'C1bar','LineStyle',linestyle,'Color',colors_vec{2});
hold on
plot(t_eval, C2, 'LineWidth', 2, 'DisplayName', 'C2','LineStyle',linestyle,'Color',colors_vec{3});
hold on
plot(t_eval, C3, 'LineWidth', 2, 'DisplayName', 'C3','LineStyle',linestyle,'Color',colors_vec{4});
hold on
plot(t_eval, C4, 'LineWidth', 2, 'DisplayName', 'C4','LineStyle',linestyle,'Color',colors_vec{5});
hold on
plot(t_eval, C5, 'LineWidth', 2, 'DisplayName', 'C5','LineStyle',linestyle,'Color',colors_vec{6});
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;

%%
subplot(rows,cols,2)
hold on
plot(t_eval, C1barC4, 'LineWidth', 2, 'DisplayName', 'C1barC4','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
ylim([0,10e-12])
%%
subplot(rows,cols,3)
hold on
plot(t_eval, C4b2, 'LineWidth', 2, 'DisplayName', 'C4b2','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, C4bC2aC3, 'LineWidth', 2, 'DisplayName', 'C4bC2aC3','LineStyle',linestyle,'Color',colors_vec{2});
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
hold off;
%%

subplot(rows,cols,4)
hold on
plot(t_eval, C4bC2a, 'LineWidth', 2, 'DisplayName', 'C4bC2a','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, C4bC2aC3bC5b, 'LineWidth', 2, 'DisplayName', 'C4bC2aC3bC5b','LineStyle',linestyle,'Color',colors_vec{2});
hold on
plot(t_eval, C4bC2aC3b, 'LineWidth', 2, 'DisplayName', 'C4bC2aC3b','LineStyle',linestyle,'Color',colors_vec{3});
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
%%
subplot(rows,cols,5)
hold on
plot(t_eval, C4a, 'LineWidth', 2, 'DisplayName', 'C4a','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, C4b, 'LineWidth', 2, 'DisplayName', 'C4b','LineStyle',linestyle,'Color',colors_vec{2});
hold on
plot(t_eval, C2b, 'LineWidth', 2, 'DisplayName', 'C2b','LineStyle',linestyle,'Color',colors_vec{3});
hold on
plot(t_eval, C3a, 'LineWidth', 2, 'DisplayName', 'C3a','LineStyle',linestyle,'Color',colors_vec{4});
hold on
plot(t_eval,C5a, 'LineWidth', 2, 'DisplayName', 'C5a','LineStyle',linestyle,'Color',colors_vec{5});
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
%%
subplot(rows,cols,6)
hold on
plot(t_eval, C6C7C8C9, 'LineWidth', 2, 'DisplayName', 'C6C7C8C9','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, MAC, 'LineWidth', 2, 'DisplayName', 'MAC','LineStyle',linestyle,'Color',colors_vec{2});
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;

%% 
figure
subplot(2,2,1)
hold on
plot(t_eval, C4, 'LineWidth', 2, 'DisplayName', 'C4','LineStyle',linestyle,'Color',colors_vec{1});
hold on
plot(t_eval, C1, 'LineWidth', 2, 'DisplayName', 'C1','LineStyle',linestyle,'Color',colors_vec{2});
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;

subplot(2,2,2)
hold on
plot(t_eval, AgAb, 'LineWidth', 2, 'DisplayName', 'AgAb','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;

subplot(2,2,3)
hold on
plot(t_eval, C4bC2aC3bC5, 'LineWidth', 2, 'DisplayName', 'C4bC2aC3bC5','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;

subplot(2,2,4)
hold on
plot(t_eval, C1bar, 'LineWidth', 2, 'DisplayName', 'C1bar','LineStyle',linestyle,'Color',colors_vec{1});
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
