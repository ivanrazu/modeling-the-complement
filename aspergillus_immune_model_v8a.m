function dydt = aspergillus_immune_model_v8a(t, y, params)
% Equations 


% Defining State Variables
Ac = y(1);
As = y(2);
Ah = y(3);
C5a = y(4);
N = y(5);
M = y(6);
D = y(7);
H = y(8);

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
Kd=params(20);

Kdn = params(21);
Kdh = params(22);
mud = params(23);
Kh = params(24);
Khh=params(25);
Khd=params(26);
muh=params(27);
Kcn=params(28);
Kcm=params(29);
KH=params(30);

Kas =params(31);
Kac = params(32);
Kahs= params(33);
Khs = params(34);
rh = params(35);
Ahmax = params(36);


% Differential Equations

% Asperg. conidia
dAc_dt =  - Ks*Ac; %- Kas*Ac/(Kac+Ac);

% Asperg. swollen conidia
dAs_dt = Kas*Ac/(Kac+Ac) - (dasm * M + dasn * N) * As - Kahs *As/(Khs+As) ; % - 0*Kah*As ;


% Asperg. hyphae
dAh_dt =  Kahs *As/(Khs+As) + rh* Ah*(1-Ah/Ahmax)  - dah * Ah/(1+1e2*Ah) * N ;

% dAh_dt =  Kahs *As/(Khs+As) + rh* Ah - 30*Ah/(1+0.1*Ah)  - dah * Ah * N ; 


% C5a anaphalytoxin
dC5a_dt = Kc * ( Kca * (As+Ah) + Kch * H) /((1 + (Kca * (As+Ah) + Kch * H)))- 0*(Kcn*N+Kcm*M)*C5a  -  muc5a*C5a;

% Neutrophilsig
dN_dt = Kn * (Knn * C5a * N + Kna * A + Knd * D)  / ( 1 + (Knn * C5a * N + Kna * A + Knd * D)  )  - mun * N;

% Macrophages
dM_dt = Km * (Kmm * C5a * M + Kma * A + Kmd * D)  / ( 1 + (Kmm * C5a * M + Kma * A + Kmd * D)  ) - mum  * M;

% Tissue damage
dD_dt = Kd*(Kdn * N + Kdh * Ah + KH*H)/( 1 + ( Kdn * N + Kdh * Ah + KH*H)) - mud*D;

% Heme
dH_dt = Kh * (Khh *Ah + Khd * D)/(1+(Khh *Ah + Khd * D)) - muh * H;

dydt = [dAc_dt,dAs_dt,dAh_dt,dC5a_dt,dN_dt, dM_dt,dD_dt, dH_dt]';

