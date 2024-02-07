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

K9 = 8 * 10^(12); % arbitrarily chosen to obtain biologically relevant trajectory
K_9 = 0.3 * 10^(-12); % arbitrarily chosen to obtain biologically relevant trajectory

KAgAb = 1 * 10^(-7);   % arbitrarily chosen to obtain biologically relevant trajectory


FH=0;
C4bp=0;
DAF=0;
CR1=0;
CR2=0;
MCP=0;

par_str={'K0', 'K1', 'K_1', 'K2', 'K3', 'K_3', 'K4', 'K5', 'K_5', 'K6', 'K7', 'K_7', 'K8', 'K9', 'K_9'...
    , 'AgAb','FH','C4bp','DAF','CR1','CR2','MCP','KAgAb'};
%%

% name_vec={'K0', 'K1', 'K_1', 'K2', 'K3', 'K_3', 'K4', 'K5', 'K_5', 'K6', 'K7', 'K_7', 'K8', 'K9', 'K_9'...
%     , 'AgAb','FH','C4bp','DAF','CR1','CR2','MCP'};
% name_vec_latex={'K0','K0', 'K1', 'K_1', 'K2', 'K3', 'K_3', 'K4', 'K5', 'K_5', 'K6', 'K7', 'K_7', 'K8', 'K9', 'K_9'...
%     , 'AgAb','FH','C4bp','DAF','CR1','CR2','MCP'};
% a_vec = [K0, K1, K_1, K2, K3, K_3, K4, K5, K_5, K6, K7, K_7, K8, K9, K_9, AgAb,FH,C4bp,DAF,CR1,CR2,MCP].*0.5;
% b_vec = [K0, K1, K_1, K2, K3, K_3, K4, K5, K_5, K6, K7, K_7, K8, K9, K_9, AgAb,FH,C4bp,DAF,CR1,CR2,MCP].*10;
%

% name_vec={'K0','K1','K_1','K2'};
% name_vec_latex={'$K_0$','$K_1$','$K_{-1}$','$K_2$'};
% 
% a_vec = [K0*0.1,K1*0.1,K_1*0.001,K2*0.00001];
% b_vec = [K0,K1*10,K_1*100,K2*10];


name_vec={'K0'};
name_vec_latex={'$K_0$'};

a_vec = [K0*0.1];
b_vec = [K0*10];

% a_vec = [K9*1e-4,K_9*5e10];
% b_vec = [K9*1e-2,K_9*1e12];

%%
% Set =1 depending on which variable want to plot 
plot_C1=0;
plot_C1bar=0;

plot_C4=1;
plot_C1barC4=0;
plot_C4a=0;
plot_C4b=0;


plot_C6C7C8C9=0;
plot_MAC=0;
%%


% Initial conditions
C1_0 = 45 * 10^(-5);
C1bar_0 = 0;
C4_0 = 310 * 10^(-5);
C1barC4_0 = 0;
C4a_0 = 0;
C4b_0 = 0;
C2_0 = 21.36 * 10^(-5);
C4b2_0 = 0;
C4bC2a_0 = 0;
C2b_0 = 0;
C3_0 = 88.9 * 10^(-5);
C4bC2aC3_0 = 0;
C4bC2ac3b_0 = 0;
C3a_0 = 0;
C5_0 = 44.44 * 10^(-5);
C4bC2aC3bC5_0 = 0;
C4bC2aC3bC5b_0 = 0;
C5a_0 = 0;
C6789_0 = 0;
MAC_0 = 0;

AgAb_0=1e2;


params = [K0, K1, K_1, K2, K3, K_3, K4, K5, K_5, K6, K7, K_7, K8, K9, K_9,FH,C4bp,DAF,CR1,CR2,MCP,KAgAb];


initial_conditions = [C1_0, C1bar_0, C4_0, C1barC4_0, C4a_0, C4b_0, C2_0, ...
    C4b2_0, C4bC2a_0, C2b_0, C3_0, C4bC2aC3_0, C4bC2ac3b_0, ...
    C3a_0, C5_0, C4bC2aC3bC5_0, C4bC2aC3bC5b_0, C5a_0, C6789_0, MAC_0,AgAb_0]/1000;

tspan = [0, 1];

% Time points for which to evaluate the solution
t_eval = linspace(tspan(1), tspan(2), 1000);

%%

np=6;

rows=1;
cols=1;


for i=1:length(a_vec)
    a=a_vec(i);
    b=b_vec(i);

    name=name_vec(i);
    name_latex=name_vec_latex(i);

    vec=par_vec(a,b,np); % vector of equally spaced values between a and b in log space

    colors=jet(length(vec));

    idx = find(ismember(par_str, name));


f=figure;
xSize = cols*10; X=xSize; ySize = rows*8;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[100 25 xSize*50 ySize*55]);
hold on;

    for j=1:length(vec)

        if idx>1 && idx<length(params)
            pars_new=[params(1:idx-1),vec(j),params(idx+1:end)];
        end

        if idx==1
            pars_new=[vec(j),params(idx+1:end)];
        end

        if idx==length(params)
            pars_new=[params(1:idx-1),vec(j)];
        end

        solv = ode23s(@classic_pathway_CS, t_eval, initial_conditions, [], pars_new);

        Ys = deval(solv, t_eval)';

        C1  = Ys(:,1);
        C1bar = Ys(:,2);
        C4 = Ys(:,3);
        C1barC4  = Ys(:,4);
        C4a=Ys(:,5);
        C4b=Ys(:,6);
        C2=Ys(:,7);
        C4b2=Ys(:,8);
        C4bC2a=Ys(:,9);
        C2b=Ys(:,10);
        C3=Ys(:,11);
        C4bC2aC3=Ys(:,12);
        C4bC2aC3b=Ys(:,13);
        C3a=Ys(:,14);
        C5=Ys(:,15);
        C4bC2aC3bC5=Ys(:,16);
        C4bC2aC3bC5b=Ys(:,17);
        C5a=Ys(:,18);
        C6C7C8C9=Ys(:,19);
        MAC=Ys(:,20);

%%%%%%%%%%%%%%%%
if plot_C1==1
        subplot(rows,cols,1)
        hold on
        set(gca,'Fontsize',30);box on;
        plot(t_eval, C1', 'Color',colors(j,:), 'LineWidth', 2);
        grid off
        hold on
        set(gca,'linewidth',2)
        xlabel('Time');
        legend('C1')
end        
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
if plot_C1bar==1
        subplot(rows,cols,2)
        hold on
        set(gca,'Fontsize',30);box on;
        plot(t_eval, C1bar', 'Color',colors(j,:), 'LineWidth', 2);
        grid off
        hold on
        set(gca,'linewidth',2)
        xlabel('Time');
        legend('C1bar')
end
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
if plot_C4==1
        subplot(rows,cols,1)
        hold on
        set(gca,'Fontsize',30);box on;
        plot(t_eval, C4', 'Color',colors(j,:), 'LineWidth', 2);
        grid off
        hold on
        set(gca,'linewidth',2)
        xlabel('Time');
        legend('C4')
end        
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
if plot_C1barC4==1
        subplot(rows,cols,4)
        hold on
        set(gca,'Fontsize',30);box on;
        plot(t_eval, C1barC4', 'Color',colors(j,:), 'LineWidth', 2);
        grid off
        hold on
        set(gca,'linewidth',2)
        xlabel('Time');
        legend('C1barC4')
end        
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
if plot_C4a==1
        subplot(rows,cols,5)
        hold on
        set(gca,'Fontsize',30);box on;
        plot(t_eval, C4a', 'Color',colors(j,:), 'LineWidth', 2);
        grid off
        hold on
        set(gca,'linewidth',2)
        xlabel('Time');
        legend('C4a')
end        
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
if plot_C4b==1
        subplot(rows,cols,6)
        hold on
        set(gca,'Fontsize',30);box on;
        plot(t_eval, C4b', 'Color',colors(j,:), 'LineWidth', 2);
        grid off
        hold on
        set(gca,'linewidth',2)
        xlabel('Time');
        legend('C4b')
end        
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
if plot_C6C7C8C9==1
        subplot(rows,cols,1)
        hold on
        set(gca,'Fontsize',30);box on;
        plot(t_eval, C6C7C8C9', 'Color',colors(j,:), 'LineWidth', 2);
        grid off
        hold on
        set(gca,'linewidth',2)
        xlabel('Time');
        legend('C6C7C8C9')
end        
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
if plot_MAC==1
        subplot(rows,cols,2)
        hold on
        set(gca,'Fontsize',30);box on;
        plot(t_eval, MAC', 'Color',colors(j,:), 'LineWidth', 2);
        grid off
        hold on
        set(gca,'linewidth',2)
        xlabel('Time');
        legend('MAC')
end

ylim([0,1e-5])



    end

    dim = [.2 .63 .9 .3];
    annotation('textbox',dim,'String',name_latex,'interpreter','latex','Fontsize',40,'FontName','Calibri' ,'Color', 'k','EdgeColor','none');


    % saveas(f,sprintf('%s%s',strcat('Senst_',char(name)),'.png'));

end






%%
function z=par_vec(A,B,np)

a=log10(A);
b=log10(B);

z=10.^(linspace(a,b,np));
end

