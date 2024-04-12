function dydt = aspergillus_immune_model_v3(t, y, params)
% Equations 


% Defining State Variables
Ac = y(1);
As = y(2);
Ah = y(3);
C5a = y(4);
N = y(5);
M = y(6);
Aif = y(7);
D = y(8);
H = y(9);

A = Ac + As + Ah ;

% Defining model parameters
dacm = params(1);      
dacn = params(2);    
Ks = params(3);     
dasm =params(4);
dasn = params(5);
Kah = params(6);
dah=params(7);
Kc =params(8);
Kca = params(9);
Kch = params(10);

muc5a = params(11);
Kn=params(12);
Knn=params(13);
Kna=params(14);
Knd =params(15);
mun =params(16);
dnc =params(17);
dns=params(18);
dnh=params(19);
Km=params(20);

Kmm=params(21);
Kma=params(22);
Kmd=params(23);
mum=params(24);
dmc =params(25);
dms=params(26);
sai=params(27);
Kai = params(28);
Kain = params(29);
Kaim = params(30);

Kaid = params(31);
muai = params(32);
Kd=params(33);
Kdn = params(34);
Kdh = params(35);
mud = params(36);
Kh = params(37);
Khh=params(38);
Khd=params(39);
muh=params(40);

Aifstar=params(41);
C5astar=params(42);

alphac=params(43);
alphas=params(44);
alphah=params(45);

% define f
F_aif = 1/(1+(Aif/Aifstar)^2) ;

% Differential Equations

% Asperg. conidia
dAc_dt =  - (dacm *  M * + dacn * N )*  Ac *  F_aif/(alphac+0*Ac) - Ks*Ac;

% Asperg. swollen conidia
dAs_dt = Ks * Ac  - (dasm * M + dasn * N) * As * F_aif/(alphas+0*As) - Kah * As;

% Asperg. hyphae
dAh_dt = Kah * As  - dah * Ah * N * F_aif/(alphah+0*0.05*Ah); 

% C5a anaphalytoxin
dC5a_dt = Kc * ( Kca * (As+Ah) + Kch * H) * F_aif/((1 + Kca * (As+Ah) + Kch * H)* F_aif)- muc5a * C5a;

% Neutrophilsig
dN_dt = Kn * (Knn * C5a * N + Kna * A + Knd * D)* F_aif  / ( 1 + (Knn * C5a * N + Kna * A + Knd * D)* F_aif  )  - mun *( 1 +  0*(dnc * Ac /( alphac+Ac) + dns *As/( alphas+As) + dnh * Ah/( alphah+Ah))  * F_aif  * 1/(1+C5a/C5astar) ) * N;

% Macrophages
dM_dt = Km * (Kmm * C5a * M + Kma * A + Kmd * D) * F_aif / ( 1 + (Kmm * C5a * M + Kma * A + Kmd * D) * F_aif ) - mum * ( 1 +  0*(dmc * Ac/(alphac+Ac) + dms * As/(alphas+As) )* F_aif ) * M;


% Ainti-inflammatory mediators
dAif_dt = sai +  Kai*(Kain *N + Kaim * M + Kaid * D)*F_aif/( 1 + (Kain *N + Kaim * M + Kaid * D)*F_aif) - muai * Aif;

% Tissue damage
dD_dt = Kd*(Kdn * N + Kdh * Ah)/( 1 + ( Kdn * N + Kdh * Ah)) - mud*D;

% Heme
dH_dt = Kh * (Khh *Ah + Khd * D)/(1+(Khh *Ah + Khd * D)) - muh * H;

dydt = [dAc_dt,dAs_dt,dAh_dt,dC5a_dt,dN_dt, dM_dt,dAif_dt,dD_dt, dH_dt]';

