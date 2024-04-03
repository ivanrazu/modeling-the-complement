function dydt = aspergillus_immune_model_v2(t, y, params)
% Equations 


% Defining State Variables
Ac = y(1);
As = y(2);
Ah = y(3);
C5a = y(4);
N = y(5);
M = y(6);
R = y(7);
Aif = y(8);
D = y(9);
H = y(10);

A = Ac + As + Ah ;

% Defining model parameters
dacm = params(1);      
dacn = params(2);    
alpha = params(3);     
Ks = params(4);     
dasm =params(5);
dasn = params(6);
beta =params(7);
Kah = params(8);
dah=params(9);
gamma=params(10);

Kc5a =params(11);
Kca = params(12);
Kch = params(13);
muc5a = params(14);
Kn=params(15);
Knn=params(16);
Kna=params(17);
Knd =params(18);
mun =params(19);
dnc =params(20);

dns=params(21);
dnh=params(22);
Km=params(23);
Kmm=params(24);
Kma=params(25);
Kmd=params(26);
mum=params(27);
dmc =params(28);
dms=params(29);
Kr = params(30);

Krn=params(31);
Krm=params(32);
mur=params(33);
drm=params(34);
drn=params(35);
sai=params(36);
Kai = params(37);
Kain = params(38);
Kaim = params(39);
Kaid = params(40);

muai = params(41);
Kd=params(42);
Kdr = params(43);
Kdn = params(44);
Kdh = params(45);
mud = params(46);
Kh = params(47);
Khh=params(48);
Khd=params(49);
muh=params(50);

Aifstar=params(51);
C5astar=params(52);

% define f
F_aif = 1/(1+(Aif/Aifstar)^2) ;

% Differential Equations

% Asperg. conidia
dAc_dt =  - (dacm *  M * + dacn * N )*  Ac * R * F_aif/(1+alpha*Ac) - Ks*Ac;

% Asperg. swollen conidia
dAs_dt = Ks * Ac  - (dasm * M + dasn * N) * As * R * F_aif/(1+beta*As) - Kah * As;

% Asperg. hyphae
dAh_dt = Kah * As  - dah * Ah * N * R * F_aif/(1+gamma*Ah); 


% C5a anaphalytoxin
dC5a_dt = Kc5a * ( Kca * (As+Ah) + Kch * H) * F_aif/((1 + Kca * A + Kch * H)* F_aif)- muc5a * C5a;

% Neutrophilsig
dN_dt = Kn * (Knn * C5a * N + Kna * A + Knd * D) * F_aif / ( 1 + (Knn * C5a * N + Kna * A + Knd * D) * F_aif )  - mun *( 1 +  (dnc * Ac /( 1+alpha*Ac) + dns *As/( 1+beta*As) + dnh * Ah/( 1+gamma*Ah) ) * R * F_aif  * 1/(1+C5a/C5astar) ) * N;

% Macrophages
dM_dt = Km * (Kmm * C5a * M + Kma * A + Kmd * D) * F_aif / ( 1 + (Kmm * C5a * M + Kma * A + Kmd * D) * F_aif ) - mum * ( 1 + (dmc * Ac/(1+alpha*Ac) + dms * As/(1+beta*As) )* R * F_aif ) * M;

% ROS
dR_dt = Kr *(Krn*N+Krm*M)/(1+(Krn*N+Krm*M)) - mur*(1 +drm*M + drn*N)*R;

% Ainti-inflammatory mediators
dAif_dt = sai +  Kai*(Kain *N + Kaim * M + Kaid * D)*F_aif/( 1 + (Kain *N + Kaim * M + Kaid * D)*F_aif) - muai * Aif;

% Tissue damage
dD_dt = Kd*(Kdr *R + Kdn * N + Kdh * Ah)/( 1 + (Kdr * R + Kdn * N + Kdh * Ah)) - mud*D;

% Heme
dH_dt = Kh * (Khh *Ah + Khd * D)/(1+(Khh *Ah + Khd * D)) - muh * H;

dydt = [dAc_dt,dAs_dt,dAh_dt,dC5a_dt,dN_dt, dM_dt,dR_dt,dAif_dt,dD_dt, dH_dt]';

