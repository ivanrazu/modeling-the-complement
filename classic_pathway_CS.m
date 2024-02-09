function dydt = classic_pathway_CS(t, y, params)
    % Equations for the activation dynamics 
    % of the initial cascade of the classic CS pathway

    % Defining State Variables
    C1 = y(1);    
    C1bar = y(2);    
    C4 = y(3);    
    C1barC4 = y(4);    
    C4a = y(5);    
    C4b = y(6);
    C2 = y(7);    
    C4b2 = y(8);    
    C4bC2a = y(9);
    C2b = y(10);    
    C3 = y(11);    
    C4bC2aC3 = y(12);
    C4bC2aC3b = y(13);
    C3a = y(14);
    C5 = y(15);    
    C4bC2aC3bC5 = y(16);
    C4bC2aC3bC5b = y(17);
    C5a = y(18);
    C6C7C8C9 = y(19);
    MAC = y(20);
    AgAb = y(21);
    
    % Defining model parameters
    K0 = params(1);
    K1 = params(2);
    K_1 = params(3);
    K2 = params(4);
    K3 = params(5);
    K_3 = params(6);
    K4 = params(7);
    K5 = params(8);
    K_5 = params(9);
    K6 = params(10);
    K7 = params(11);
    K_7 = params(12);
    K8 = params(13);
    K9 = params(14);
    K_9 = params(15);

    FH=params(16);
    C4bp=params(17);
    DAF=params(18);
    CR1=params(19);
    CR2=params(20);
    MCP=params(21);

    KAgAb=params(22);
    KC1inh=params(23);
    C1inh=params(24);


    % Differential Equations
    dC1_dt = - K0 * C1 * AgAb;

    dC1bar_dt = K0 * C1 * AgAb - K1 * C1bar * C4 + K_1 * C1barC4 + K2 * C1barC4 - KC1inh * C1bar * C1inh;

    dC4_dt = -K1 * C1bar * C4 + K_1 * C1barC4;

    dC1barC4_dt = K1 * C1bar * C4 - K_1 * C1barC4 - K2 * C1barC4;
    
    dC4a_dt = K2 * C1barC4;
    
    dC4b_dt = -K3 * C4b * C2 + K_3 * C4b2 + K2 * C1barC4;
    
    dC2_dt = -K3 * C4b * C2 + K_3 * C4b2;

    dC4b2_dt = K3 * C4b * C2 - K_3 * C4b2 - K4 * C4b2;

    dC4bC2a_dt = K4 * C4b2 - K5 * C4bC2a * C3 + K_5 * C4bC2aC3 - FH - C4bp - DAF;

    dC2b_dt = K4 * C4b2;
    
    dC3_dt = -K5 * C4bC2a * C3 + K_5 * C4bC2aC3;

    dC4bC2aC3_dt = K5 * C4bC2a * C3 - K_5 * C4bC2aC3 - K6 * C4bC2aC3 + FH + C4bp;
    
    dC4bC2aC3b_dt = K6 * C4bC2aC3 - K7 * C4bC2aC3b * C5 + K_7 * C4bC2aC3bC5 - FH- C4bp - DAF - CR1 - CR2 -MCP;
    
    dC3a_dt = K6 * C4bC2aC3;

    dC5_dt = -K7 * C4bC2aC3b * C5 + K_7 * C4bC2aC3bC5;
    
    dC4bC2aC3bC5_dt = K7 * C4bC2aC3b * C5 - K_7 * C4bC2aC3bC5 - K8 * C4bC2aC3bC5 + FH;
    
    dC4bC2aC3bC5b_dt = K8 * C4bC2aC3bC5 - K9 * C4bC2aC3bC5b * C6C7C8C9 + K_9 * MAC;
    
    C5a_dt = K8 * C4bC2aC3bC5;
    
    dC6C7C8C9_dt = -K9 * C4bC2aC3bC5b * C6C7C8C9 + K_9 * MAC;
    
    dMAC_dt = K9 * C4bC2aC3bC5b * C6C7C8C9 - K_9 * MAC;

    dAgAb_dt = -KAgAb*AgAb;
    
    dydt = [dC1_dt, dC1bar_dt,dC4_dt,dC1barC4_dt,dC4a_dt,dC4b_dt,dC2_dt, ...
            dC4b2_dt,dC4bC2a_dt,dC2b_dt,dC3_dt,dC4bC2aC3_dt,dC4bC2aC3b_dt,...
             dC3a_dt, dC5_dt, dC4bC2aC3bC5_dt,dC4bC2aC3bC5b_dt,C5a_dt,dC6C7C8C9_dt,dMAC_dt,dAgAb_dt]';

