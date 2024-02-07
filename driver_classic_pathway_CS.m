clear

% Define parameters
K0 = 0.5 * 10^(1); % arbitrarily chosen to obtain biologically relevant trajectory

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
% 
% AgAb = 2.5 * 10^7;   % arbitrarily chosen to obtain biologically relevant trajectory

KAgAb = 1 * 10^(-7);

FH=0;
C4bp=0;
DAF=0;
CR1=0;
CR2=0;
MCP=0;

% K0=K0*0.000000001;


% K9=K9*1e-4;
% K_9=K_9*5e10;


params = [K0, K1, K_1, K2, K3, K_3, K4, K5, K_5, K6, K7, K_7, K8, K9, K_9,FH,C4bp,DAF,CR1,CR2,MCP,KAgAb];

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
t_span = [0, 1.5];

% Call solver
sol = ode23s(@classic_pathway_CS, t_span, initial_conditions, [], params);

% Time points for which to evaluate the solution
t_eval = linspace(t_span(1), t_span(2), 1000);

% Evaluate solution
sol_values = deval(sol, t_eval);

% Plot the results
cols=3;
rows=2;

figure;

subplot(rows,cols,1)
hold on
plot(t_eval, sol_values(1, :), 'LineWidth', 2, 'DisplayName', 'C1');
hold on;
plot(t_eval, sol_values(2, :), 'LineWidth', 2, 'DisplayName', 'C1bar');
hold on
plot(t_eval, sol_values(7, :), 'LineWidth', 2, 'DisplayName', 'C2');
hold on
plot(t_eval, sol_values(11, :), 'LineWidth', 2, 'DisplayName', 'C3');
hold on
plot(t_eval, sol_values(3, :), 'LineWidth', 2, 'DisplayName', 'C4');
hold on
plot(t_eval, sol_values(15, :), 'LineWidth', 2, 'DisplayName', 'C5');
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
hold off;

%%
% % Plot the results
% figure;

subplot(rows,cols,2)
plot(t_eval, sol_values(13, :), 'LineWidth', 2, 'DisplayName', 'C4bC2aC3b');
hold on
plot(t_eval, sol_values(4, :), 'LineWidth', 2, 'DisplayName', 'C1barC4');
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
hold off;
%%
subplot(rows,cols,3)
plot(t_eval, sol_values(8, :), 'LineWidth', 2, 'DisplayName', 'C4b2');
hold on
plot(t_eval, sol_values(12, :), 'LineWidth', 2, 'DisplayName', 'C4bC2aC3');
hold on
plot(t_eval, sol_values(16, :), 'LineWidth', 2, 'DisplayName', 'C4bC2aC3bC5');
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
hold off;
%%

% figure;
subplot(rows,cols,4)
plot(t_eval, sol_values(9, :), 'LineWidth', 2, 'DisplayName', 'C4bC2a');
hold on
plot(t_eval, sol_values(17, :), 'LineWidth', 2, 'DisplayName', 'C4bC2aC3bC5b');
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
hold off;
%%
% figure;
subplot(rows,cols,5)
plot(t_eval, sol_values(5, :), 'LineWidth', 2, 'DisplayName', 'C4a');
hold on
plot(t_eval, sol_values(6, :), 'LineWidth', 2, 'DisplayName', 'C4b');
hold on
plot(t_eval, sol_values(10, :), 'LineWidth', 2, 'DisplayName', 'C2b');
hold on
plot(t_eval, sol_values(14, :), 'LineWidth', 2, 'DisplayName', 'C3a');
hold on
plot(t_eval, sol_values(18, :), 'LineWidth', 2, 'DisplayName', 'C5a');
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
hold off;
%%
% figure;
subplot(rows,cols,6)
hold on
plot(t_eval, sol_values(19, :), 'LineWidth', 2, 'DisplayName', 'C6C7C8C9');
hold on
plot(t_eval, sol_values(20, :), 'LineWidth', 2, 'DisplayName', 'MAC');
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
hold off;

%% 
figure
plot(t_eval, sol_values(21, :), 'LineWidth', 2, 'DisplayName', 'AgAb');
xlabel('Time');
ylabel('Population');
legend('Location', 'Best');
title('Classic Pathway Dynamics');
grid on;
hold off;