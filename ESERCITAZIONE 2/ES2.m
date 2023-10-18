%% ESERCITAZIONE 2

clear
close all
clc


V2=5;             %[l]
k01=1.2;          %[h^-1]
k02=1.2;          %[h^-1]
k21_lin=2.2;      %[h^-1]
D_es=500;         %[mg] bolo per via orale
Vmax=110;         %[mg/h]
Km=50;            %[mg]
t=0:0.07:10;      %[h]

%% A

% IMPOSTO L'ODE PER L'ANDAMENTO NON LINEARE

[T_nl,Q_nl] = ode45(@(t,q) ASSORBIMENTO_NO_LIN_MM(q, k01, k02, Km, Vmax),[0 10],[D_es 0]); 

% IMPOSTO L'ODE PER L'ANDAMENTO LINEARE

[T_l,Q_l] = ode45(@(t,q) ASSORBIMENTO_LIN(q, k01, k02, k21_lin),[0 10],[D_es 0]); 

figure(1)

subplot(2,1,1)
plot(T_nl,Q_nl(:,1),'-o',T_nl,Q_nl(:,2),'-o'), grid on
title('ANDAMENTO QUANTITA Michaelis-Menten')
legend('COMP 1','COMP 2')
xlabel('[h]') 
ylabel('[mg]') 

subplot(2,1,2)
plot(T_l,Q_l(:,1),'-o',T_l,Q_l(:,2),'-o'), grid on
title('ANDAMENTO QUANTITA CINETICA LINEARE')
legend('COMP 1','COMP 2')
xlabel('[h]') 
ylabel('[mg]') 

k21_MM = (Vmax)./((Km)+Q_nl(:,1));

%%% MICHAELIS MENTEN

figure(2)

subplot(3,1,1)

plot(T_nl,k21_MM)
title('VELOCITA MM')
xlabel('tempo [h]')
ylabel('velocita [1/h]')

% FLUSSO MM

subplot(3,1,2)
plot(T_nl,Q_nl(:,1).*k21_MM)
title('FLUSSO')
xlabel('tempo [ore]')
ylabel('flusso [mg/h]')

% ANDAMENTO FLUSSO IN FUZNIONE DI Q

subplot(3,1,3)
plot(Q_nl(:,1),Q_nl(:,1).*k21_MM)
title('FLUSSO formato MM')
xlabel('quantità [mg]')
ylabel('flusso [mg/h]')


%%% LINEARE
figure(3)

subplot(3,1,1)

yline(k21_lin)
title('VELOCITA LINEARE')
xlabel('tempo [h]')
ylabel('velocita [1/h]')

% FLUSSO MM

subplot(3,1,2)
plot(T_l,Q_l(:,1).*k21_lin)
title('FLUSSO')
xlabel('tempo [ore]')
ylabel('flusso [mg/h]')

% ANDAMENTO FLUSSO IN FUZNIONE DI Q

subplot(3,1,3)
plot(Q_l(:,1),Q_l(:,1).*k21_lin)
title('FLUSSO LINEARE')
xlabel('quantità [mg]')
ylabel('flusso [mg/h]')


% Che differenze ci sono tra gli andamenti di questo sistema e quelli del 
% modello lineare corrispondente ottenuto ponendo k21 = 2.2 ore−1?


% NEL CASO NON LINEARE LA CINETICA DEL FARMACO E' PIU' RAPIDA RISPETTO AL
% CASO LINEARE. QUESTO SI PUO' OSSERVARE DALL'ANDAMENTO DELLE QUANTITA' DI
% FARMACO CHE DI FATTO SONO PIU' ALTE NEL CASO LINEARE.

%% B) 
% Cosa varia se Vmax raddopia?
% e se km si riduce della meta?

% PER VMAX CHE RADDOPPIA MI ASPETTO UN ASSORBIMENTO PIU' RAPIDO MENTRE PER
% km DIMEZZATO MI ASPETTO SEMPRE UN ASSORBIMENTO PIU' RAPIDO

% RADDOPPIO Vmax  [mg/h]

[T_nl_B1,Q_nl_B1] = ode45(@(t,q) ASSORBIMENTO_NO_LIN_MM(q, k01, k02, Km, Vmax*2),[0 10],[D_es 0]); 

% DIMEZZO km [mg]

[T_nl_B2,Q_nl_B2] = ode45(@(t,q) ASSORBIMENTO_NO_LIN_MM(q, k01, k02, Km/2, Vmax),[0 10],[D_es 0]); 

figure(4)

subplot(3,1,1)
plot(T_nl,Q_nl(:,1),'-o',T_nl,Q_nl(:,2),'-o'), grid on
title('ANDAMENTO Michaelis-Menten')
legend('COMP 1','COMP 2')
xlabel('[h]') 
ylabel('[mg]') 

subplot(3,1,2)
plot(T_nl_B1,Q_nl_B1(:,1),'-o',T_nl_B1,Q_nl_B1(:,2),'-o'), grid on
title('ANDAMENTO Michaelis-Menten Vmax*2')
legend('COMP 1','COMP 2')
xlabel('[h]') 
ylabel('[mg]') 

subplot(3,1,3)
plot(T_nl_B2,Q_nl_B2(:,1),'-o',T_nl_B2,Q_nl_B2(:,2),'-o'), grid on
title('ANDAMENTO Michaelis-Menten km/2')
legend('COMP 1','COMP 2')
xlabel('[h]') 
ylabel('[mg]') 

figure(41)

k21_MM_1 = (2*Vmax)./((Km)+Q_nl(:,1));

subplot(3,1,1)

plot(T_nl,k21_MM)
title('VELOCITA MM V RADDOPPIATO')
xlabel('tempo [h]')
ylabel('velocita [1/h]')
hold on 
plot(T_nl,k21_MM_1)


% FLUSSO MM

subplot(3,1,2)
plot(T_nl,Q_nl(:,1).*k21_MM)
title('FLUSSO V RADDOPPIATO')
xlabel('tempo [ore]')
ylabel('flusso [mg/h]')
hold on
plot(T_nl,Q_nl(:,1).*k21_MM_1)


% ANDAMENTO FLUSSO IN FUZNIONE DI Q

subplot(3,1,3)
plot(Q_nl(:,1),Q_nl(:,1).*k21_MM)
title('FLUSSO formato MM V RADDOPPIATO')
xlabel('quantità [mg]')
ylabel('flusso [mg/h]')
hold on
plot(Q_nl(:,1),Q_nl(:,1).*k21_MM_1)


figure(42)

k21_MM_2 = (Vmax)./((Km/2)+Q_nl(:,1));

subplot(3,1,1)

plot(T_nl,k21_MM)
title('VELOCITA MM Km DIMEZZATO')
xlabel('tempo [h]')
ylabel('velocita [1/h]')
hold on 
plot(T_nl,k21_MM_2)


% FLUSSO MM

subplot(3,1,2)
plot(T_nl,Q_nl(:,1).*k21_MM)
title('FLUSSO Km DIMEZZATO')
xlabel('tempo [ore]')
ylabel('flusso [mg/h]')
hold on
plot(T_nl,Q_nl(:,1).*k21_MM_2)


% ANDAMENTO FLUSSO IN FUZNIONE DI Q

subplot(3,1,3)
plot(Q_nl(:,1),Q_nl(:,1).*k21_MM)
title('FLUSSO formato MM Km DIMEZZATO')
xlabel('quantità [mg]')
ylabel('flusso [mg/h]')
hold on
plot(Q_nl(:,1),Q_nl(:,1).*k21_MM_2)

%% C)

% Cosa cambia riducendo a 50 mg la dose somministrata?
% Si confronti questo sistema con il sistema lineare corrispondente 
% ottenuto ponendo k21 = 2.2 ore−1?

D_esC = 50; %[mg]

[T_nl_C1,Q_nl_C1] = ode45(@(t,q) ASSORBIMENTO_NO_LIN_MM(q, k01, k02, Km, Vmax),[0 10],[D_esC 0]); 

[T_l_C1,Q_l_C1] = ode45(@(t,q) ASSORBIMENTO_LIN(q, k01, k02, k21_lin),[0 10],[D_esC 0]); 

figure(5)
subplot(2,2,1)
plot(T_nl_C1,Q_nl_C1(:,1),'-o',T_nl_C1,Q_nl_C1(:,2),'-o'), grid on
title('ANDAMENTO QUANTITA Michaelis-Menten 50mg')
legend('COMP 1','COMP 2')
xlabel('[h]') 
ylabel('[mg]') 

subplot(2,2,3)
plot(T_l_C1,Q_l_C1(:,1),'-o',T_l_C1,Q_l_C1(:,2),'-o'), grid on
title('ANDAMENTO QUANTITA CINETICA LINEARE k21=2.2')
legend('COMP 1','COMP 2')
xlabel('[h]') 
ylabel('[mg]') 

subplot(2,2,2)
plot(T_nl,Q_nl(:,1),'-o',T_nl,Q_nl(:,2),'-o'), grid on
title('ANDAMENTO QUANTITA Michaelis-Menten 50mg')
legend('COMP 1','COMP 2')
xlabel('[h]') 
ylabel('[mg]') 

subplot(2,2,4)
plot(T_l,Q_l(:,1),'-o',T_l,Q_l(:,2),'-o'), grid on
title('ANDAMENTO QUANTITA CINETICA LINEARE k21=2.2')
legend('COMP 1','COMP 2')
xlabel('[h]') 
ylabel('[mg]') 


figure(52)

k21_MM_C = (Vmax)./((Km)+Q_nl_C1(:,1));

subplot(3,1,1)

plot(T_nl_C1,k21_MM_C)
title('VELOCITA MM E LINEARE DOSE = 50')
xlabel('tempo [h]')
ylabel('velocita [1/h]')
hold on 
yline(k21_lin,'r')


% FLUSSO MM

subplot(3,1,2)
plot(T_nl_C1,Q_nl_C1(:,1).*k21_MM_C)
title('FLUSSO MM E LIN DOSE = 50')
xlabel('tempo [ore]')
ylabel('flusso [mg/h]')
hold on
plot(T_l_C1,Q_l_C1(:,1)*k21_lin)


% ANDAMENTO FLUSSO IN FUZNIONE DI Q

subplot(3,1,3)
plot(Q_nl_C1(:,1),Q_nl_C1(:,1).*k21_MM_C)
title('FLUSSO formato MM E LINEARE DOSE = 50')
xlabel('quantità [mg]')
ylabel('flusso [mg/h]')
hold on
plot(Q_l_C1(:,1),Q_l_C1(:,1)*k21_lin)


%% D)

Q1=1;
Q2=5;
Q3=20;
[T_nl_1,Q_nl_1] = ode45(@(t,q) ASSORBIMENTO_NO_LIN_HILL(q, k01, k02, Km, Vmax,Q1),[0 10],[D_es 0]); 
[T_nl_2,Q_nl_2] = ode45(@(t,q) ASSORBIMENTO_NO_LIN_HILL(q, k01, k02, Km, Vmax,Q2),[0 10],[D_es 0]); 
[T_nl_3,Q_nl_3] = ode45(@(t,q) ASSORBIMENTO_NO_LIN_HILL(q, k01, k02, Km, Vmax,Q3),[0 10],[D_es 0]); 


figure(6)
subplot(3,1,1)
k21_hill_1=(Vmax.*(Q_nl_1(:,1).^(Q1-1)))./((Km).^(Q1)+Q_nl_1(:,1).^(Q1));
plot(Q_nl_1(:,1),k21_hill_1.*Q_nl_1(:,1),'-o'), grid on
title('ANDAMENTO QUANTITA DI FARMACO HILL Q = 1 (MENTEN)')
legend('COMP 1')
ylabel('FLUSSO [mg/h]') 
xlabel('QUANTITA[mg]') 

subplot(3,1,2)
k21_hill_2=(Vmax.*(Q_nl_2(:,1).^(Q2-1)))./((Km).^(Q2)+Q_nl_2(:,1).^(Q2));
plot(Q_nl_2(:,1),k21_hill_2.*Q_nl_2(:,1),'-o'), grid on
title('ANDAMENTO QUANTITA DI FARMACO HILL Q = 5')
legend('COMP 1')
ylabel('FLUSSO [mg/h]') 
xlabel('QUANTITA[mg]') 

subplot(3,1,3)
k21_hill_3=(Vmax.*(Q_nl_3(:,1).^(Q3-1)))./((Km).^(Q3)+Q_nl_3(:,1).^(Q3));
plot(Q_nl_3(:,1),k21_hill_3.*Q_nl_3(:,1),'-o'), grid on
title('ANDAMENTO QUANTITA DI FARMACO HILL Q = 20')
legend('COMP 1')
ylabel('FLUSSO [mg/h]') 
xlabel('QUANTITA[mg]') 

%% E)

% %%%%%%%%%%%%%%%%%%%%%% BIODISPONIBILITA' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Biodisponibilità MENTEN 

F_MM=k21_MM./(k21_MM+k01);
ka_menten=k21_MM+k01;

% Biodisponibilità HILL CON Q=5
Q_b=5;
k21_hill=(Vmax.*(Q_nl_1(:,1).^(Q_b-1)))./((Km).^(Q_b)+Q_nl_1(:,1).^(Q_b));

F2=k21_hill./(k21_hill+k01);
ka_hill=k21_hill+k01;

k21=2.2;

F_lin=k21./(k21+k01);
ka_lineare=k21+k01;

%% confronto DOSI RIPETUTE LINEARE e NON LINEARE

% caso LINEARE con dose orale cfr. es 5 esercitazione 1

A1 = [(-k01-k21_lin) 0;k21_lin -k02];
Bor = [1;0];
C1 = [0 1/V2];
D1 = [0];

SYSRL = ss(A1,Bor,C1,D1);

% con RL = ripetute lineare
AI_RL = 500;          % [mg]
DI_R = 4;             % [ore] [ogni quanto gli do la dose]
NI_RL = 12;           % [scelto da me]


tRL = linspace(0,DI_R,40);
xRL = [0 0];
YRL = [];
XRL = [];
TRL = [];
% questo ciclo mi restituisce la concentrazione e mi va bene
% perchè è quello che chiede il testo
for i = 1:NI_RL
   [ylRL,tlRL,xlRL] = initial(SYSRL,xRL,tRL);
   [yfRL,tfRL,xfRL] = impulse(SYSRL*AI_RL,tRL);
   YRL = [YRL; ylRL+yfRL];
   XRL = [XRL; xlRL+xfRL];
   TRL = [TRL tRL+DI_R*(i-1)];
   xRL = xfRL(end,:) + xlRL(end,:);
end



% caso NON LINEARE con dosi ripetute
% considero una non linearità di tipo M-M
% con RNL = ripetute non lineare
dose_R = 500;     % [mg]
DI_R = 2;         % [ore] [ogni quanto gli do la dose]
NI_RNL = 12;      % [scelto da me]

tRNL = linspace(0,DI_R,40);

x1RNL = dose_R;     % cond iniziali compartimento 1
x2RNL = 0;          % cond iniziali compartimento 2
TRNL_MM = [];
QRNL_MM = [];

for j = 1:NI_RNL
    [T1RNL,Q1RNL] = ode45(@(t,q) ASSORBIMENTO_NO_LIN_MM...
        (q, k01, k02, Km, Vmax),tRNL,[x1RNL x2RNL]); 

    TRNL_MM = [TRNL_MM;T1RNL];
    QRNL_MM = [QRNL_MM;Q1RNL];
    x1RNL = Q1RNL(1, end) + dose_R;    % compartimento 1
    x2RNL = Q1RNL(2, end);             % compartimento 2

    tRNL = linspace(j*DI_R,(j+1)*DI_R,40);
end 

% questo ciclo mi restituisce la quantità 
% devo allora calcolare la concentrazione perchè è quello
% che mi chiede il testo C=Q/V
CRNL_MM = QRNL_MM/V2;



% caso NON LINEARE con dosi ripetute
% considero una non linearità di tipo HILL
x3RNL = dose_R;     % cond iniziali compartimento 1
x4RNL = 0;          % cond iniziali compartimento 2
TRNL_HILL = [];
QRNL_HILL = [];
QR = 10;           % parametro della hill
                    % scelgo un Q "basso" e diverso da 1 al fine di 
                    % sottolineare eventuali differenze col caso M-M

tRNL2 = linspace(0,DI_R,40);
          

for k = 1:NI_RNL
   [T2RNL,Q2RNL] = ode45(@(t,q) ASSORBIMENTO_NO_LIN_HILL...
       (q, k01, k02, Km, Vmax, QR),tRNL2,[x3RNL x4RNL]); 

    TRNL_HILL = [TRNL_HILL;T2RNL];
    QRNL_HILL = [QRNL_HILL;Q2RNL];
    x3RNL = Q2RNL(1, end) + dose_R;    % compartimento 1
    x4RNL = Q2RNL(2, end);             % compartimento 2

    tRNL2 = linspace(k*DI_R,(k+1)*DI_R,40);
end 

% questo ciclo mi restituisce la quantità 
% devo allora calcolare la concentrazione perchè è quello
% che mi chiede il testo C=Q/V
CRNL_HILL = QRNL_HILL/V2;



% confronto le risposte LIN - MM - HILL
% aka plotto i grafici assieme
figure(10)
subplot(3,1,1)
plot(TRL,YRL)
grid on
title('LINEARE')
xlabel('tempo [ore]')
ylabel('concentrazione [mg/l]')
subplot(3,1,2)
plot(TRNL_MM,CRNL_MM)
grid on
title('NON lineare M-M')
xlabel('tempo [ore]')
ylabel('concentrazione [mg/l]')
subplot(3,1,3)
plot(TRNL_HILL,CRNL_HILL)
grid on
title('NON lineare HILL [Q=10]')
xlabel('tempo [ore]')
ylabel('concentrazione [mg/l]')






