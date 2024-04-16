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
  
Ks = params(1);     
dasm =params(2);
dasn = params(3);
Kah = params(4);
dah=params(5);
Kc =params(6);
Kca = params(7);
Kch = params(8);

muc5a = params(9);
Kn=params(10);
Knn=params(11);
Kna=params(12);
Knd =params(13);
mun =params(14);
Km=params(15);
Kmm=params(16);
Kma=params(17);
Kmd=params(18);

mum=params(19);
sai=params(20);
Kai = params(21);
Kain = params(22);
Kaim = params(23);
Kaid = params(24);
muai = params(25);
Kd=params(26);
Kdn = params(27);
Kdh = params(28);

mud = params(29);
Kh = params(30);
Khh=params(31);
Khd=params(32);
muh=params(33);
Aifstar=params(34);

mu=params(35);
sigma=params(36);
Kcn=params(37);
Kcm=params(38);

% define f
F_aif = 1/(1+(Aif/Aifstar)^2) ;

% Differential Equations

% Asperg. conidia
% dAc_dt =  - (dacm *  M * + dacn * N )*  Ac *  F_aif- Ks*1/sqrt(2*pi*sigma^2)*exp(-(t-mu)^2/(2*sigma^2));
dAc_dt =  - Ks*1/sqrt(2*pi*sigma^2)*exp(-(t-mu)^2/(2*sigma^2));

% Asperg. swollen conidia
dAs_dt = Ks*1/sqrt(2*pi*sigma^2)*exp(-(t-mu)^2/(2*sigma^2))  - (dasm * M + dasn * N) * As * F_aif -As/sqrt(2*pi*sigma^2)*exp(-(t-(mu+4))^2/(2*sigma^2)) ;

% Asperg. hyphae
dAh_dt =  As/sqrt(2*pi*sigma^2)*exp(-(t-(mu+4))^2/(2*sigma^2)) + Kah *As  - dah * Ah * N * F_aif; 

% C5a anaphalytoxin
dC5a_dt = Kc * ( Kca * (As+Ah) + Kch * H)* F_aif /((1 + (Kca * (As+Ah) + Kch * H))* F_aif)- muc5a*(1+Kcn*N+Kcm*M) * C5a;

% Neutrophilsig
dN_dt = Kn * (Knn * C5a * N + Kna * A + Knd * D)* F_aif  / ( 1 + (Knn * C5a * N + Kna * A + Knd * D)* F_aif  )  - mun * N;

% Macrophages
dM_dt = Km * (Kmm * C5a * M + Kma * A + Kmd * D) * F_aif / ( 1 + (Kmm * C5a * M + Kma * A + Kmd * D) * F_aif ) - mum  * M;


% Ainti-inflammatory mediators
dAif_dt = sai +  Kai*(Kain *N + Kaim * M + Kaid * D)*F_aif/( 1 + (Kain *N + Kaim * M + Kaid * D)*F_aif) - muai * Aif;

% Tissue damage
dD_dt = Kd*(Kdn * N + Kdh * Ah)/( 1 + ( Kdn * N + Kdh * Ah)) - mud*D;

% Heme
dH_dt = Kh * (Khh *Ah + Khd * D)/(1+(Khh *Ah + Khd * D)) - muh * H;

dydt = [dAc_dt,dAs_dt,dAh_dt,dC5a_dt,dN_dt, dM_dt,dAif_dt,dD_dt, dH_dt]';

