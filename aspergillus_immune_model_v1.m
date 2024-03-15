function dydt = aspergillus_immune_model_v1(t, y, params)
% Equations 


% Defining State Variables
Ac = y(1);
Ah = y(2);
C5a = y(3);
N = y(4);
M = y(5);
R = y(6);
Aif = y(7);
D = y(8);

A = Ac + Ah ;

% Defining model parameters
Kac = params(1);      
Acmax = params(2);    
dacm = params(3);     
dacn = params(4);     
Kah =params(5);
Ahmax = params(6);
dahn =params(7);
Kc5a = params(8);
Aifstar=params(9);
muc5a=params(10);
Kmm=params(11);
mum=params(12);
Knn=params(13);
mun =params(14);
Krn =params(15);
Krm =params(16);
mur=params(17);
sai=params(18);
Kain=params(19);
Kaim=params(20);
muai=params(21);
Kdr=params(22);
Kdn=params(23);
mud =params(24);
dn=params(25);
dm = params(26);
Kai=params(27);
drm=params(28);
drn=params(29);
Kd=params(30);
alpha=params(31);
beta=params(32);
Kna = params(33);
Kma = params(34);
Kn = params(35);
Km = params(36);
Kr = params(37);

% Differential Equations

% Asperg. conidia
dAc_dt = Kac * Ac * (1-Ac/Acmax) - 0.1*Ac/(1 + Ac)   - 1/(1+(Aif/Aifstar)^2)*dacm  * Ac * M * R/(1+alpha*Ac) - 1/(1+(Aif/Aifstar)^2)*dacn * Ac * N * R/(1+alpha*Ac);

% Asperg. hyphae
dAh_dt = Kah * Ac  - dahn  * Ah * N * R/(1+beta*Ah) ;

% C5a anaphalytoxin
dC5a_dt = Kc5a * A/(1+Aif/Aifstar)^2 - muc5a * C5a;

% Neutrophils
dN_dt = Kn*(Knn* C5a * N + Kna * A)/(1+(Aif/Aifstar)^2) - mun * (1+dn*Ac*R )*N ;

% Macrophages
dM_dt = Km*(Kmm * C5a * M + Kma * A)/(1+(Aif/Aifstar)^2) - mum * (1+dm*Ac*R)*M;

% ROS
dR_dt = Kr *(Krn*N+Krm*M)/(1+(Krn*N+Krm*M)) - mur*(1 +drm*M + drn*N)*R;

% Ainti-inflammatory mediators
dAif_dt = sai +  Kai*(Kain *N + Kaim * M)/( 1 + (Kain *N + Kaim * M)) - muai * Aif;

% Tissue damage
dD_dt = Kd*(Kdr *R + Kdn * N )/( 1 + (Kdr * R + Kdn * N )) - mud*D;




dydt = [dAc_dt,dAh_dt,dC5a_dt,dN_dt, dM_dt,dR_dt,dAif_dt,dD_dt]';

