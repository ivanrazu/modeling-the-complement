function dydt = Complement_System_mechanistic_v2(t, y, params)
% Equations for the activation dynamics
% of the initial cascade of the classic CS pathway

% includes c3b, disticntion in C3 convertases between CP,LP and AP
% includes Asperg conidia and hyphae

% Defining State Variables
C1 = y(1);
C2 = y(2);
C3 = y(3);
C4 = y(4);
C5 = y(5);
C3convCPLP = y(6);
C3convAP = y(7);
C5conv = y(8);
C3a = y(9);
C5a = y(10);
C3b = y(11);
C3H2O = y(12); % C3 hydrolized
Ac = y(13);   % Aspergillus conidia
Ah = y(14);   % Aspergillus hyphae
Am = y(15); % Alveolar Macrophages
N = y(16);  % Neutrophils
C1inh =y(17) ;
FB = y(18); % Factor B
FH = y(19); % Factor H
FI = y(20); % Factor I
P = y(21);  % Properdin


C3conv = C3convCPLP + C3convAP;
A = Ac + Ah ;

% Defining model parameters

Kc1 = params(1);           Kc2 = params(2);          Kc3 = params(3);
Kc4 = params(4);           Kc5 = params(5);          sc1 = params(6);
sc2 = params(7);           sc3 = params(8);          sc4 = params(9);
sc5 = params(10);          muc1 = params(11);        muc2 = params(12);
muc3 = params(13);         muc4 = params(14);        muc5 = params(15);
C1star = params(16);       C2star = params(17);      C3star = params(18);
C4star = params(19);       C5star = params(20);      Kc3convcpcl = params(21);
Kcp = params(22);          Klp = params(23);         C1inhstar = params(24);
muc3convcplp = params(25); Kc3convap = params(26);   C3H2Ostar = params(27);
muc3convap = params(28);   Kc5conv = params(29);     Kc5convhs = params(30);
muc5conv = params(31);     Kc3acplp = params(32);    Kc3aap = params(33);
muc3a = params(34);        Kc5a = params(35);        muc5a = params(36);
Kc3bcplp = params(37);     Kc3bh2ofb = params(38);   Kc3bap = params(39);
FHstar = params(40);       FIstar = params(41);      muc3b = params(42);
Kc3h2o = params(43);       muc3h2o = params(44);     KAc = params(45);
Acmax = params(46);        dacm = params(47);        dacn = params(48);
KAh = params(49);          Ahmax = params(50);       dahm = params(51);
dahn = params(52);         Kam = params(53);         Kamc3a = params(54);
Kamc5a = params(55);       Amhs = params(56);        muam = params(57);
dam = params(58);          Kn = params(59);          Knc3a = params(60);
Knc5a = params(61);        Nhs = params(62);         mun = params(63);
                                                     dan = params(64);

sc1inh = params(65);
muc1inh = params(66);
sfb = params(67);
mufb = params(68);
sfh = params(69);
mufh = params(70);
sfi = params(71);
mufi = params(72);
Kp =params(73);
mup = params(74);
Pstar =params(75);

% Differential Equations

dC1_dt = sc1 + Kc1 * C1 * A/(1+C1/C1star) - muc1 * C1;
dC2_dt = sc2 + Kc2 * C2 * A/(1+C2/C2star) - muc2 * C2;
dC3_dt = sc3 + Kc3 * C3 * A/(1+C3/C3star) - muc3 * C3;
dC4_dt = sc4 + Kc4 * C4 * A/(1+C4/C4star) - muc4 * C4;
dC5_dt = sc5 + Kc5 * C5 * A/(1+C5/C5star) - muc5 * C5;

% C3 convertase from Classical Pathway (CP) and Lectin Pathway (LP)
dC3convCPLP_dt = Kc3convcpcl * A *( ( Kcp * C1 * C2 * C4) + (Klp * C2 *  C4) )/( 1 + C1inh/C1inhstar ) -  muc3convcplp * C3convCPLP;

% C3 convertase frrom Alternative Pathway (AP)
% dC3convAP_dt = Kc3convap * A *( C3H2O/(1+C3H2O/C3H2Ostar) ) -  muc3convap * C3convAP;
dC3convAP_dt = Kc3convap  * ( C3H2O * FB )/(1+FH/FHstar) -  muc3convap * C3convAP/(1+P/Pstar)*(1+FH);

% C5 convertase
dC5conv_dt = Kc5conv * ( C3b * C3conv)/( Kc5convhs + ( C3b * C3conv))  -  muc5conv * C5conv;

dC3a_dt = Kc3acplp *  C3 * C3convCPLP + Kc3aap * C3 * C3convAP  - muc3a * C3a;

dC5a_dt = Kc5a *  C5 * C5conv  - muc5a * C5a;

% C3b can be rpoduced from CP-LP an frrom AP
dC3b_dt = Kc3bcplp * C3 * C3convCPLP + (Kc3bh2ofb * C3H2O * FB  + Kc3bap* C3convAP * C3) / ( (1+FH/FHstar)*(1+FI/FIstar)  ) - muc3b * C3b * (1+FH);

dC3H2O_dt = Kc3h2o * C3 - muc3h2o * C3H2O;



% Asperg. conidia
% if Ac<0
%     Ac=0;
%     dAc_dt=0;
% else
    dAc_dt = KAc * Ac * (1-Ac/Acmax) -dacm  * Ac * Am * C3b - dacn * Ac * N * C3b;
% end


% Asperg. hyphae
% if Ah<0
%     Ah=0;
%     dAh_dt=0;
% else
    dAh_dt = KAh * Ac * (1-Ac/Ahmax) - dahm  * Ah * Am * C3b - dahn * Ah * N * C3b;
% end

% Alveolar macrophages
if Am<0
    Am=0;
    dAm_dt=0;
else
    dAm_dt = Kam *A*(Kamc3a * C3a + Kamc5a * C5a )/(Amhs + (Kamc3a * C3a + Kamc5a * C5a )) - muam * (1 + dam * A * C3b)*Am ;
end


% Neutrophils
if N<0
    N=0;
    dN_dt=0;
else
    dN_dt = Kn * A *(Knc3a * C3a + Knc5a * C5a)/(Nhs+(Knc3a * C3a + Knc5a * C5a)) - mun *(1 + dan * A * C3b) *N ;
end

% Inhibitors
dC1inh_dt = sc1inh - muc1inh * C1inh;
dFB_dt = sfb - mufb * FB;
dFH_dt = sfh - mufh * FH;
dFI_dt = sfi - mufi* FI;
dP_dt= Kp - mup * P;

dydt = [dC1_dt,dC2_dt,dC3_dt,dC4_dt,dC5_dt,...
    dC3convCPLP_dt,dC3convAP_dt,dC5conv_dt,...
    dC3a_dt,dC5a_dt,dC3b_dt,dC3H2O_dt,dAc_dt,...
    dAh_dt,dAm_dt, dN_dt, dC1inh_dt, dFB_dt, dFH_dt,dFI_dt,dP_dt]';

