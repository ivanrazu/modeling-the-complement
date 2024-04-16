function dydt = aspergillus_immune_model_v6(t, y, params)
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
Km=params(17);
Kmm=params(18);
Kma=params(19);
Kmd=params(20);

mum=params(21);
sai=params(22);
Kai = params(23);
Kain = params(24);
Kaim = params(25);
Kaid = params(26);
muai = params(27);
Kd=params(28);
Kdn = params(29);
Kdh = params(30);

mud = params(31);
Kh = params(32);
Khh=params(33);
Khd=params(34);
muh=params(35);
Aifstar=params(36);

mu=3;
sigma=1;

% define f
F_aif = 1/(1+(Aif/Aifstar)^2) ;

% Differential Equations

% Asperg. conidia
dAc_dt =  - (dacm *  M * + dacn * N )*  Ac *  F_aif- Ks*1/sqrt(2*pi*sigma^2)*exp(-(t-mu)^2/(2*sigma^2));

% Asperg. swollen conidia
dAs_dt = Ks*1/sqrt(2*pi*sigma^2)*exp(-(t-mu)^2/(2*sigma^2))  - (dasm * M + dasn * N) * As * F_aif - As/sqrt(2*pi*sigma^2)*exp(-(t-(mu+4))^2/(2*sigma^2)) ;

% Asperg. hyphae
dAh_dt =  As/sqrt(2*pi*sigma^2)*exp(-(t-(mu+4))^2/(2*sigma^2)) + Kah *As  - dah * Ah * N * F_aif; 

% C5a anaphalytoxin
% dC5a_dt = Kc * ( Kca * (As+Ah) + Kch * H) * F_aif/((1 + Kca * (As+Ah) + Kch * H)* F_aif)- muc5a*(1+0.05*N+0.03*M) * C5a;
dC5a_dt = Kc * ( Kca * (As+Ah) + Kch * H)* F_aif /((1 + (Kca * (As+Ah) + Kch * H))* F_aif)- muc5a*(1+0.8*N+0.8*M) * C5a;

% Neutrophilsig
dN_dt = Kn * (Knn * C5a * N + Kna * A + Knd * D)* F_aif  / ( 500 + (Knn * C5a * N + Kna * A + Knd * D)* F_aif  )  - mun * N;

% Macrophages
dM_dt = Km * (Kmm * C5a * M + Kma * A + Kmd * D) * F_aif / ( 100 + (Kmm * C5a * M + Kma * A + Kmd * D) * F_aif ) - mum  * M;


% Ainti-inflammatory mediators
dAif_dt = sai +  Kai*(Kain *N + Kaim * M + Kaid * D)*F_aif/( 1 + (Kain *N + Kaim * M + Kaid * D)*F_aif) - muai * Aif;

% Tissue damage
dD_dt = Kd*(Kdn * N + Kdh * Ah)/( 1 + ( Kdn * N + Kdh * Ah)) - mud*D;

% Heme
dH_dt = Kh * (Khh *Ah + Khd * D)/(1+(Khh *Ah + Khd * D)) - muh * H;

dydt = [dAc_dt,dAs_dt,dAh_dt,dC5a_dt,dN_dt, dM_dt,dAif_dt,dD_dt, dH_dt]';

