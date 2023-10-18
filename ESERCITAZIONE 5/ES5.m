%% ESERCITAZIONE 5 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% ES 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clc
clear
close all
% CASO E >> S
% la funzione 'ENZIMA' modellizza la CINETICA ENZIMATICA a DUE STADI 

% UNITA' DI MISURA
% 1 Millimolare [mM] = 0,001 Mole per litro [mol/l]


% DATI

k1e = 1; % k1 velocita di formazione del composto [mM^-1*sec^-1]
k2e = 1; % k2 velocita di dissociazione del composto e formazione del prodotto [sec^-1]
k_1e = 1; % k_1 velocita di dissociazione del composto [sec^-1]
k_2e = 1; % k_2 velocita di ricombinazione prodotto enzima [mM^-1*sec^-1]

% CONDIZIONI INIZIALI 

s0 = 1;
e0 = 100;
p0 = 0;
c0 = 0;
t = [0:0.001:0.15];

[tout,cout] = ode45( @(t,x)ENZIMA(x ,k1e, k2e, k_1e, k_2e ), t, [s0,e0,c0,p0]) ;

% PRODOTTO allo SS e TEMPO in cui lo raggiungo
% MASSIMA CONCENTRAZIONE DEL PRODOTTO E TEMPO DI RAGGIUNGIMENTO CON FIND

p_ss = max(cout(:,4));
t_p_ss = (tout(76)); 


% GRAFICO ANDAMENTI DELLE CONCENTRAZIONI

figure('Name','ANDAMENTO CONCENTRAZIONI E>>S')

subplot(2,2,1)
plot(tout,cout(:,1))
title('Substrato')   % DA 1 passa a 0.01
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]') 
 
subplot(2,2,2)
plot(tout,cout(:,2))
title('Enzima')       % DA 100 passa a 99.02
xlabel('Tempo [sec]') 
ylabel('Concentrazione [mM]') 

subplot(2,2,3)
plot(tout,cout(:,3))
title('Complesso')   % DA 0 passa a 0.97
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')

subplot(2,2,4)
plot(tout,cout(:,4))
axis([0 0.15 0 0.015]) 
hold on
yline(p_ss,'r') % CONCENTRAZIONE PRODOTTO ALLO SS
hold on
xline(t_p_ss,'r') % TEMPO RAGGIUNGIMENTO CONCENTRAZIONE ALLO SS
title('Prodotto')  % DA 0 passa a 0.01
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% andamento grafico CRESCITA PRODOTTO e DECRESCITA SUBSTRATO
% per vedere se effettivamente è un andamento exp, lo plotto su scala semilog
% P vedo che ho un andamento crescente ma non posso dire che siano exp
% S vedo che ho un andamento exp

figure('Name','scala SEMILOGARITMICA')

subplot(2,1,1)
semilogy(tout,cout(:,4))
title('Prodotto')
xlabel('Tempo [sec]')
ylabel('Concentrazione semilog')

subplot(2,1,2)
semilogy(tout,cout(:,1))
title('Substrato')
xlabel('Tempo [sec]')
ylabel('Concentrazione semilog')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% VELOCITA DI FORMAZIONE DEL PRODOTTO

% vp = k2*E-k_2*P*E
vp = k2e.*cout(:,3)-k_2e.*cout(:,2).*cout(:,4);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in QUANTO TEMPO si raggiunge lo STEADY STATE
% allo STEADY STATE la velocità di formazione del prodotto si annulla
% considero il secondo degli indici perchè il primo è la condizione
% iniziale di 'partenza' della velocità

ind_ss = find(vp<0.001);
t_ss = (t(ind_ss(2)));
figure('Name','Vp con E>>S')
plot(t,vp)
hold on
xline(t_ss,'r')
title('Velocita di formazione del prodotto E>>S')
xlabel('Tempo [sec]')
ylabel('Velocita [mM/sec]')
legend('Vp','tempo SS')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% ES 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% CASO S >> E

% dati che variano
% nb. i k rimangono invariati

e02 = 1;
s02 = 100;
t_se = [0:0.1:400]; 

[tout2, cout2] = ode45( @(t,x)ENZIMA(x ,k1e, k2e, k_1e, k_2e ), t_se, [s02,e02,c0,p0]) ;


% PRODOTTO allo SS e TEMPO in cui lo raggiungo
% MASSIMA CONCENTRAZIONE DEL PRODOTTO E TEMPO DI RAGGIUNGIMENTO CON FIND


p_ss2 = max(cout2(:,4));
t_p_ss2 = (tout2(end));


% GRAFICO ANDAMENTI DELLE CONCENTRAZIONI

figure('Name','Andamenti delle cocentrazioni S>>E')

subplot(2,2,1)
plot(tout2, cout2(:,1))
title('Substrato')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')
 
subplot(2,2,2)          % istantaneamente l'enzima si lega 
plot(tout2, cout2(:,2)) % e non è più enzima libero 
axis([-0.2 1 0 1])      % ma va a far parte di c
title('Enzima')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')

subplot(2,2,3)
plot(tout2, cout2(:,3)) 
axis([-0.2 1 0 1]) 
title('Complesso')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')

subplot(2,2,4)
plot(tout2, cout2(:,4))
axis([0 450 0 60]) 
hold on
yline(p_ss2,'g')
hold on
xline(t_p_ss2,'g')
title('Prodotto')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% SEMILOG

figure('Name','scala SEMILOGARITMICA')

subplot(2,1,1)
semilogy(tout2,cout2(:,4))
title('Prodotto')
xlabel('Tempo [sec]')
ylabel('Concentrazione semilog')

subplot(2,1,2)
semilogy(tout2,cout2(:,1))
title('Substrato')
xlabel('Tempo [sec]')
ylabel('Concentrazione semilog')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% VELOCITA DI FORMAZIONE DEL PRODOTTO

% vp = k2*E-k_2*P*E
vp2 = k2e.*cout2(:,3)-k_2e.*cout2(:,2).*cout2(:,4);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in QUANTO TEMPO si raggiunge lo STEADY STATE

ind_ss2 = find(vp2<0.001);
t_ss2 = (t_se(ind_ss2(2)));

figure('Name','Vp con S>>E')

plot(t_se,vp2)
axis([-50 400 0 1.1]) 
hold on
xline(t_ss2,'r')
title('Velocita di formazione del prodotto')
xlabel('Tempo [sec]')
ylabel('Velocita [mM/sec]')
legend('Vp','tempo SS')

%% %%%%%%%%%%%%%%%%%%%%%%%%% ES 1 VS ES 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ANDAMENTI CONCENTRAZIONI

figure('Name','CONFRONTO Andamenti delle concentrazioni')

subplot(4,2,1)                    % S1
plot(tout,cout(:,1))
title('Substrato E>>S')
xlabel('Tempo [sec]')
ylabel('Conce [mM]')
 
subplot(4,2,3)                    % E1
plot(tout,cout(:,2))
title('Enzima E>>S')
xlabel('Tempo [sec]')
ylabel('Conc [mM]')

subplot(4,2,5)                    % C1
plot(tout,cout(:,3))
title('Complesso E>>S')
xlabel('Tempo [sec]')
ylabel('Conc [mM]')

subplot(4,2,7)                    % P1
plot(tout,cout(:,4))
axis([0 0.15 0 0.015]) 
hold on
yline(p_ss,'g')
hold on
xline(t_p_ss,'g')
title('Prodotto E>>S')
xlabel('Tempo [sec]')
ylabel('Conc [mM]')

subplot(4,2,2)                    % S2
plot(tout2, cout2(:,1))
title('Substrato S>>E')
xlabel('Tempo [sec]')
ylabel('Conc [mM]')
 
subplot(4,2,4)                    % E2
plot(tout2, cout2(:,2))
axis([-0.2 1 0 1]) 
title('Enzima S>>E')
xlabel('Tempo [sec]')
ylabel('Conc [mM]')

subplot(4,2,6)                    % C2
plot(tout2, cout2(:,3)) 
axis([-0.2 1 0 1]) 
title('Complesso S>>E')
xlabel('Tempo [sec]')
ylabel('Conc [mM]')

subplot(4,2,8)                    % P2
plot(tout2, cout2(:,4))
axis([0 450 0 60]) 
hold on
yline(p_ss2,'g')
hold on
xline(t_p_ss2,'g')
title('Prodotto S>>E')
xlabel('Tempo [sec]')
ylabel('Conc [mM]')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% CONFRONTO LE VELOCITA DI FORMAZIONE DEL PRODOTTO

% NELLA PRIMA FASE LA VELOCITA 1 E' MOLTO PIU' VELOCE DELLA 2
% E' POSSIBILE FARE APPROSSIMAZIONE ALL'EQUILIBRIO E STATO QUASI
% STAZIONARIO

ind_ss = find(vp<0.001);
t_ss = (t(ind_ss(2)));
ind_ss2 = find(vp2<0.15);
t_ss2 = (t_se(ind_ss2(2)));


figure('Name','CONFRONTO VP')

subplot(2,1,1)                    % VP1
plot(t,vp)
hold on
xline(t_ss,'r')
title('E >> S')
xlabel('Tempo [sec]')
ylabel('Velocita [mM/sec]')

subplot(2,1,2)                    % VP2
plot(t_se,vp2)
axis([-50 400 0 1.1])
hold on
xline(t_ss2,'r')
title('S >> E')
xlabel('Tempo [sec]')
ylabel('Velocita [mM/sec]')


%% %%%%%%%%%%%%%%%%%%%% ES 3 carbonic-anhydrase %%%%%%%%%%%%%%%%%%%%%%%%%%%


% SPESSO LA RICOMBINAZIONE E' ASSENTE COME NEL CASO DELL'IDRATAZIONE DELLA
% CO2 DURANTE IL TRASFERIMENTO DEI TESSUTI AL SANGUE E DA QUESTO ALL'ARIA
% ALVEOLARE CATALIZZATO DALL'ENZIMA CARBONIC-ANHYDRASE

% CASO S>>E

% DATI

k1co2 = 75e3;
k2co2 = 600e3;
k_1co2 = 75;
k_2co2 = 0; % !!!!!!!!! IRREVERSIBILE FORMAZIONE DEL PRODOTTO !!!!!!!!!!!
e0co2 = 1;
s0co2 = 100;
p0co2 = 0;
c0co2 = 0;
t_co2 = 0:0.00001:5e-4; 

[toutco2, coutco2] = ode45( @(t,x)ENZIMA(x ,k1co2, k2co2, k_1co2, k_2co2 ), ...
    t_co2, [s0co2,e0co2,c0co2,p0co2]) ;


% PRODOTTO allo SS e TEMPO in cui lo raggiungo
% MASSIMA CONCENTRAZIONE DEL PRODOTTO E TEMPO DI RAGGIUNGIMENTO CON FIND


p_ssco2 = max(coutco2(:,4));
t_p_ssco2 = (toutco2(38));


% GRAFICO ANDAMENTI DELLE CONCENTRAZIONI


figure('Name','Andamenti delle cocentrazioni CO2')

subplot(2,2,1)                    % S CO 
plot(toutco2, coutco2(:,1))
title('Substrato')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')


subplot(2,2,2)                    % E CO 
plot(toutco2, coutco2(:,2)) 
title('Enzima')                   % LIBERO A REGIME
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')

subplot(2,2,3)                    % C CO 
plot(toutco2, coutco2(:,3))
title('Complesso')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')

subplot(2,2,4)                    % P CO 
plot(toutco2, coutco2(:,4))
hold on
yline(p_ssco2,'g')
hold on
xline(t_p_ssco2,'g')
title('Prodotto')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% andamento grafico CRESCITA PRODOTTO e DECRESCITA SUBSTRATO
% per vedere se effettivamente è un andamento exp, lo plotto su scala semilog
% nel grafico di crescita del prodotto potrebbe esserci una prima parte di
% andamento exp (da 0 a 10)

% nel grafico di decrescita del substrato vedo una chiara divisione in
% rette decrescenti quindi presuppongo un andamento exp

% DUE COMPONENTI

figure('Name','scala SEMILOGARITMICA')

subplot(2,1,1)
semilogy(toutco2,coutco2(:,4))
title('Prodotto')
xlabel('Tempo [sec]')
ylabel('Concentrazione semilog')

subplot(2,1,2)
semilogy(toutco2,coutco2(:,1))
title('Substrato')
xlabel('Tempo [sec]')
ylabel('Concentrazione semilog')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% VELOCITA' DI FORMAZIONE DEL PRODOTTO

% vp = k2*C-k_2*P*E
vpco2 = k2co2.*coutco2(:,3)-k_2co2.*coutco2(:,2).*coutco2(:,4);

ind_ssco2= find(vpco2<0.1);
t_ssco2 = (t_co2(ind_ssco2(2)));

figure('Name','Vp CO2')
plot(t_co2,vpco2)
hold on
xline(t_ssco2,'r')
title('Vp CO2')
xlabel('Tempo [sec]')
ylabel('Velocita [mM/sec]')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% D) E E) 


% SCALA TEMPORALE
% si evidenzia bene sia dal grafico dell'andamento di p
% sia da quello di velocità di formazione del prodotto che l'ordine di
% grandezza a cui va a regime è 10^-4 secondi
% si notano due scale temporali: una veloce e una lenta


% tau0 concentrazione COMPLESSO
% difficile comparare su scale temporali cosi piccole, ma risultano dello
% stesso ordine di grandezza quindi direi di si; l'errore potrebbe essere
% dato dal fatto che l'asse dei tempi è discretizzato dalla ode45
% RMB. tau=1/K
% t = 5*tau
% velocità di formazione del complesso
% non scrivo la parte con k_2 perchè tanto è = 0

vcco2 = (k1co2*e0co2.*coutco2(:,1)) -((k_1co2+k2co2+k1co2.*coutco2(:,1)).*coutco2(:,3)); 

% tau e tempo che stimo con la formula
Tau0 = 1/((k1co2*s0co2)+k_1co2+k2co2);

ttt=1/(k1co2*s0co2);
TTT=(s0co2+(k_1co2/k1co2))/(k2co2*e0co2);
tTau0 = 5*Tau0;

% tau e tempo che ricavo dai dati
ind_tau0 = find( vcco2 < 0.01*max(vcco2), 1 );
t_tau0 = toutco2(ind_tau0);
Tau_0 = t_tau0/5;
figure('Name','Andamenti della concentrazione del COMPLESSO CO2')

plot(t_co2,vcco2)
xline(tTau0,'m')
hold on 
xline(t_tau0,'k')
hold on 
xline(TTT,'r')
hold on 
xline(ttt,'g')
hold on 
xline(Tau0,'b')

title('Complesso')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Michaelis Menten e Lineweaver-Burk Plot

% come secondo me ha più senso
% RMB. M-M LWBP
% asse x = s0 asse x = 1/s0
% asse y = Vp asse y = 1/Vp

% IMPOSTO I PARAMETRI DELLA MM ------> QUASI ALLO STEADY STATE

Kmco2 = (k_1co2+k2co2)/(k1co2);

Vmaxco2 = k2co2*e0co2;

S0CO2 = 0:1:1000;

Vmmco2 = (Vmaxco2.*S0CO2)./(Kmco2+S0CO2);

figure('Name','Michaelis-Menten')

subplot (2,1,1)
plot(S0CO2,Vmmco2)
axis([-100 1000 0 1e6])
hold on 
yline(Vmaxco2,'r')
hold on
xline(Kmco2,'m')
title('M-M co2')
xlabel('s0')
ylabel('Vp [mM/sec]')
legend('curva','Vmax','Km')

subplot(2,1,2)
plot(1./S0CO2,1./Vmmco2)
title('Lineweaver-Burk Plot co2')
xlabel('1/s0')
ylabel('1/Vp')


%% %%%%%%%%%%%%%%%%%% ES 4 INIBIZIONE ENZIMATICA %%%%%%%%%%%%%%%%%%%%%%%%%%

% simulazione dei diversi tipi di inibizione al variare di s0
% ATTENZIONE: uso la ode23s perchè il tutor sostiene che sia più stabile

% DATI

s0in = 100;        %[mM]
e0in = 1;          %[mM]
c0in = 0;
p0in = 0;
c20in = 0;
p20in = 0;
k1in = 75e3;
k2in = 600e3;
k_1in = 75;
k_2in = 0;
k3in = 75e2; 
k_3in = 75; 
k4in = 600e3; 
k_4in = 0; 

% il tempo lo scelgo io opportunamente a posteriori 

Tend = 0.004; 
t = 0:0.0001:Tend;

% IMPOSTO LE ODE NEI 3 DIVERSI CASI

[Tno_compbase, Pno_compbase] = ode23s(@(t,p) NON_COMPETITIVA... 
    (p, k1in, k_1in, k2in, k_2in, k3in, k_3in, k4in, k_4in),...
    [0 Tend],[s0in e0in c0in p0in 100 c20in p20in]);

[Tcompbase, Pcompbase] = ode23s(@(t,p) COMPETITIVA...
    (p, k1in, k_1in, k2in, k3in, k_3in),...
    [0 Tend],[s0in e0in c0in p0in 100 c20in]);

[Tincompbase, Pincompbase] = ode23s(@(t,p) INCOMPETITIVA...
    (p, k1in, k_1in, k2in, k3in, k_3in), ...
    [0 Tend],[s0in e0in c0in p0in 100 p20in]);

figure('Name','3 casi di inibizione (base)')


subplot(2,4,1)                                   % S 
plot(Tno_compbase,Pno_compbase(:,1),'c')
hold on
plot(Tcompbase,Pcompbase(:,1),'b')
hold on
plot(Tincompbase,Pincompbase(:,1),'g')
title('substrato [S]')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')
hold on


subplot(2,4,3)                                   % ES 
plot(Tno_compbase,Pno_compbase(:,3),'c')
hold on
plot(Tcompbase,Pcompbase(:,3),'b')
hold on
plot(Tincompbase,Pincompbase(:,3),'g')
title('complesso [ES]')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')
hold on


subplot(2,4,4)                                   % P 
plot(Tno_compbase,Pno_compbase(:,4),'c')
hold on
plot(Tcompbase,Pcompbase(:,4),'b')
hold on
plot(Tincompbase,Pincompbase(:,4),'g')
title('prodotto [P]')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')
hold on


subplot(2,4,5)                                   % I
plot(Tno_compbase,Pno_compbase(:,5),'c')
hold on
plot(Tcompbase,Pcompbase(:,5),'b')
hold on
plot(Tincompbase,Pincompbase(:,5),'g')
title('inibitore [I]')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')
hold on


subplot(2,4,6)                                   % EI
plot(Tno_compbase,Pno_compbase(:,6),'c')
hold on
plot(Tcompbase,Pcompbase(:,6),'b')
title('complesso 2 [EI]')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')
hold on



subplot(2,4,7)                                   % EIS 
plot(Tno_compbase,Pno_compbase(:,7),'c')
hold on
plot(Tincompbase,Pincompbase(:,6),'g')
title('prodotto 2 [EIS]')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')


subplot(2,4,2)                                   % E 
plot(Tno_compbase,Pno_compbase(:,2),'c')
hold on
plot(Tcompbase,Pcompbase(:,2),'b')
hold on
plot(Tincompbase,Pincompbase(:,2),'g')
title('enzima [E]')
xlabel('Tempo [sec]')
ylabel('Concentrazione [mM]')
hold on
legend({'non competitiva','competitiva','incompetitiva'}, ...
 'Position',[0.747,0.109,0.156,0.344]) 
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%%%%%%%%%%%%%%%%%%% INIBIZONE NON COMPETITIVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCELGO DIVERSI VALORI PER L'INIBITORE

i0in = [1 10 20 50 100 200 500 1000];

sNO_COMP = cell(1,length(i0in));
eNO_COMP = cell(1,length(i0in));
cNO_COMP = cell(1,length(i0in));
pNO_COMP = cell(1,length(i0in));
iNO_COMP = cell(1,length(i0in));
eiNO_COMP = cell(1,length(i0in));
eisNO_COMP = cell(1,length(i0in));


for i = 1:length(i0in)
 [Tno_comp, Pno_comp] = ode23s(@(t,p) NON_COMPETITIVA...
     (p, k1in, k_1in, k2in, k_2in, k3in, k_3in, k4in, k_4in), ...
     [0 Tend],[s0in e0in c0in p0in i0in(i) c20in p20in]);
 
 
 sNO_COMP{1,i} = [Tno_comp,Pno_comp(:,1)];
 eNO_COMP{1,i} = [Tno_comp,Pno_comp(:,2)];
 cNO_COMP{1,i} = [Tno_comp,Pno_comp(:,3)];
 pNO_COMP{1,i} = [Tno_comp,Pno_comp(:,4)];
 iNO_COMP{1,i} = [Tno_comp,Pno_comp(:,5)];
 eiNO_COMP{1,i} = [Tno_comp,Pno_comp(:,6)];
 eisNO_COMP{1,i} = [Tno_comp,Pno_comp(:,7)];
end  

% AD OGNI CICLO VIENE CREATA UNA POSIZIONE NEI CELL ARRAY CHE CONTIENE 
% IL VETTORE DI OUTPUT DELLA ODE

figure('Name','Non Competitiva')

hold on

for i = 1:length(i0in)
 

 subplot(2,4,1)                                   % S 
 plot(sNO_COMP{1,i}(:,1),sNO_COMP{1,i}(:,2))
 title('substrato [S]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,2)                                   % E 
 plot(eNO_COMP{1,i}(:,1),eNO_COMP{1,i}(:,2))
 title('enzima [E]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on
 
 subplot(2,4,3)                                   % ES 
 plot(cNO_COMP{1,i}(:,1),cNO_COMP{1,i}(:,2))
 title('complesso [ES]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,4)                                   % P 
 plot(pNO_COMP{1,i}(:,1),pNO_COMP{1,i}(:,2))
 title('prodotto [P]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,5)                                   % I 
 plot(iNO_COMP{1,i}(:,1),iNO_COMP{1,i}(:,2))
 title('inibitore [I]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,6)                                   % EI 
 plot(eiNO_COMP{1,i}(:,1),eiNO_COMP{1,i}(:,2))
 title('complesso 2 [EI]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on
 
 subplot(2,4,7)                                   % EIS 
 plot(eisNO_COMP{1,i}(:,1),eisNO_COMP{1,i}(:,2))
 title('prodotto 2 [EIS]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on
 

 hold on
 legend({'i0=1','i0=10','i0=20','i0=50','i0=100','i0=200', ...
 'i0=500','i0=1000'},'Position',[0.747,0.109,0.156,0.344]) 

end 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%%%%%%%%%%%%%%%%%%%%%%% INIBIZONE COMPETITIVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sCOMP = cell(1,length(i0in));
eCOMP = cell(1,length(i0in));
cCOMP = cell(1,length(i0in));
pCOMP = cell(1,length(i0in));
iCOMP = cell(1,length(i0in));
eiCOMP = cell(1,length(i0in));

for i = 1:length(i0in)
 [Tcomp, Pcomp] = ode23s(@(t,p) COMPETITIVA(p, k1in, k_1in, k2in, k3in, k_3in), ...
 [0 Tend],[s0in e0in c0in p0in i0in(i) c20in]);
 
 
 sCOMP{1,i} = [Tcomp,Pcomp(:,1)];
 eCOMP{1,i} = [Tcomp,Pcomp(:,2)];
 cCOMP{1,i} = [Tcomp,Pcomp(:,3)];
 pCOMP{1,i} = [Tcomp,Pcomp(:,4)];
 iCOMP{1,i} = [Tcomp,Pcomp(:,5)];
 eiCOMP{1,i} = [Tcomp,Pcomp(:,6)];
 
end 

figure('Name','Competitiva')

hold on

for i = 1:length(i0in)
 
 subplot(2,4,1)                                   % S
 plot(sCOMP{1,i}(:,1),sCOMP{1,i}(:,2))
 title('substrato [S]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,2)                                   % E
 plot(eCOMP{1,i}(:,1),eCOMP{1,i}(:,2))
 title('enzima [E]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,3)                                   % ES
 plot(cCOMP{1,i}(:,1),cCOMP{1,i}(:,2))
 title('complesso [ES]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,4)                                   % P
 plot(pCOMP{1,i}(:,1),pCOMP{1,i}(:,2))
 title('prodotto [P]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,5)                                   % I
 plot(iCOMP{1,i}(:,1),iCOMP{1,i}(:,2))
 title('inibitore [I]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,6)                                   % EI
 plot(eiCOMP{1,i}(:,1),eiCOMP{1,i}(:,2))
 title('complesso 2 [EI]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 hold on
 legend({'i0=1','i0=10','i0=20','i0=50','i0=100','i0=200', ...
 'i0=500','i0=1000'},'Position',[0.747,0.109,0.156,0.344]) 

end 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%%%%%%%%%%%%%%%%%%%%% INIBIZONE INCOMPETITIVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sINCOMP = cell(1,length(i0in));
eINCOMP = cell(1,length(i0in));
cINCOMP = cell(1,length(i0in));
pINCOMP = cell(1,length(i0in));
iINCOMP = cell(1,length(i0in));
eisINCOMP = cell(1,length(i0in));

for i = 1:length(i0in)
 [Tinco, Pinco] = ode23s(@(t,p) INCOMPETITIVA...
     (p, k1in, k_1in, k2in, k3in, k_3in), ...
 [0 Tend],[s0in e0in c0in p0in i0in(i) p20in]);
 
 sINCOMP{1,i} = [Tinco,Pinco(:,1)];
 eINCOMP{1,i} = [Tinco,Pinco(:,2)];
 cINCOMP{1,i} = [Tinco,Pinco(:,3)];
 pINCOMP{1,i} = [Tinco,Pinco(:,4)];
 iINCOMP{1,i} = [Tinco,Pinco(:,5)];
 eisINCOMP{1,i} = [Tinco,Pinco(:,6)];
 
end 
figure('Name','Incompetitiva')

hold on
for i = 1:length(i0in)
 

 subplot(2,4,1)                                   % S
 plot(sINCOMP{1,i}(:,1),sINCOMP{1,i}(:,2))
 title('substrato [S]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,2)                                   % E
 plot(eINCOMP{1,i}(:,1),eINCOMP{1,i}(:,2))
 title('enzima [E]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,3)                                   % ES
 plot(cINCOMP{1,i}(:,1),cINCOMP{1,i}(:,2))
 title('complesso [ES]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,4)                                   % P
 plot(pINCOMP{1,i}(:,1),pINCOMP{1,i}(:,2))
 title('prodotto [P]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,5)                                   % I
 plot(iINCOMP{1,i}(:,1),iINCOMP{1,i}(:,2))
 title('inibitore [I]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,6)                                   % EIS
 plot(eisINCOMP{1,i}(:,1),eisINCOMP{1,i}(:,2))
 title('prodotto 2 [EIS]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 hold on
 legend({'i0=1','i0=10','i0=20','i0=50','i0=100','i0=200', ...
 'i0=500','i0=1000'},'Position',[0.747,0.109,0.156,0.344]) 

 
end 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% CONFRONTO LE VELOCITA DI FORMAZIONE DEL PRODOTTO

% al variare di i0 (considero i 3 diversi casi)

VpCOMP = cell(1,length(i0in));
VpNO_COMP = cell(1,length(i0in));
VpINCOMP = cell(1,length(i0in));

for i = 1:length(i0in)
 % COMPETITIVA
 VpCOMP{1,i} = k2in.*cCOMP{1,i}; 
 
 % NON COMPETITIVA
 VpNO_COMP{1,i} = k2in*cNO_COMP{1,i};
 
 % INCOMPETITIVA
 VpINCOMP{1,i} = k2in*cINCOMP{1,i}; 
end

figure('Name','Confronto Vp nei 3 casi')

hold on

for i = 1:length(i0in)
 

 subplot(3,2,1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPETITIVA
 plot(cCOMP{1,i}(:,1),VpCOMP{1,i}(:,2))
 title('competitiva')
 xlabel('Tempo [sec]')
 ylabel('Velocità [mM*sec-1]')
 hold on
 
 subplot(3,2,3) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NON COMPETITIVA
 plot(cNO_COMP{1,i}(:,1),VpNO_COMP{1,i}(:,2))
 title('non competitiva')
 xlabel('Tempo [sec]')
 ylabel('Velocità [mM*sec-1]')
 hold on

 subplot(3,2,5) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INCOMPETITIVA
 plot(cINCOMP{1,i}(:,1),VpINCOMP{1,i}(:,2))
 title('incompetitiva')
 xlabel('Tempo [sec]')
 ylabel('Velocità [mM*sec-1]')
 hold on
 
 legend({'i0=1','i0=10','i0=20','i0=50','i0=100','i0=200', ...
 'i0=500','i0=1000'},'Position',[0.747,0.109,0.156,0.344]) 

 
end 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALISI DI SENSITIVITA'

% come cambia il tutto al variare di k3
% considero i 3 differenti casi

k3new = [25e2 50e2 75e2 150e2 250e2 500e2];
s0new = 100; 
i0new = 100;


sNO_COMP_K = cell(1,length(k3new));
eNO_COMP_K = cell(1,length(k3new));
cNO_COMP_K = cell(1,length(k3new));
pNO_COMP_K = cell(1,length(k3new));
iNO_COMP_K = cell(1,length(k3new));
eiNO_COMP_K = cell(1,length(k3new));
eisNO_COMP_K = cell(1,length(k3new));


sCOMP_K = cell(1,length(k3new));
eCOMP_K = cell(1,length(k3new));
cCOMP_K = cell(1,length(k3new));
pCOMP_K = cell(1,length(k3new));
iCOMP_K = cell(1,length(k3new));
eiCOMP_K = cell(1,length(k3new));

sINCOMP_K = cell(1,length(k3new));
eINCOMP_K = cell(1,length(k3new));
cINCOMP_K = cell(1,length(k3new));
pINCOMP_K = cell(1,length(k3new));
iINCOMP_K = cell(1,length(k3new));
eisINCOMP_K = cell(1,length(k3new));

for i = 1:length(k3new)
 %%%%%%%%%%%%%%%%%%%% NON COMPETITIVA analisi di SENSITIVITA'

 [Tinn, Pinn] = ode23s(@(t,p) NON_COMPETITIVA...
     (p, k1in, k_1in, k2in, k_2in, k3new(i), k_3in, k4in, k_4in), ...
     [0 Tend],[s0new e0in c0in p0in i0new c20in p20in]);
 

 sNO_COMP_K{1,i} = [Tinn,Pinn(:,1)];
 eNO_COMP_K{1,i} = [Tinn,Pinn(:,2)];
 cNO_COMP_K{1,i} = [Tinn,Pinn(:,3)];
 pNO_COMP_K{1,i} = [Tinn,Pinn(:,4)];
 iNO_COMP_K{1,i} = [Tinn,Pinn(:,5)];
 eiNO_COMP_K{1,i} = [Tinn,Pinn(:,6)];
 eisNO_COMP_K{1,i} = [Tinn,Pinn(:,7)];


 %%%%%%%%%%%%%%%%%%%% COMPETITIVA analisi di SENSITIVITA'

 [Tcompn, Pcompn] = ode23s(@(t,p) COMPETITIVA...
     (p, k1in, k_1in, k2in, k3new(i), k_3in), ...
 [0 Tend],[s0new e0in c0in p0in i0new c20in]);
 

 sCOMP_K{1,i} = [Tcompn,Pcompn(:,1)];
 eCOMP_K{1,i} = [Tcompn,Pcompn(:,2)];
 cCOMP_K{1,i} = [Tcompn,Pcompn(:,3)];
 pCOMP_K{1,i} = [Tcompn,Pcompn(:,4)];
 iCOMP_K{1,i} = [Tcompn,Pcompn(:,5)];
 eiCOMP_K{1,i} = [Tcompn,Pcompn(:,6)];


 %%%%%%%%%%%%%%%%%%%% INCOMPETITIVA analisi di SENSITIVITA'
 [Tincon, Pincon] = ode23s(@(t,p) INCOMPETITIVA...
     (p, k1in, k_1in, k2in, k3new(i), k_3in), ...
 [0 Tend],[s0new e0in c0in p0in i0new p20in]);
 

 sINCOMP_K{1,i} = [Tincon,Pincon(:,1)];
 eINCOMP_K{1,i} = [Tincon,Pincon(:,2)];
 cINCOMP_K{1,i} = [Tincon,Pincon(:,3)];
 pINCOMP_K{1,i} = [Tincon,Pincon(:,4)];
 iINCOMP_K{1,i} = [Tincon,Pincon(:,5)];
 eisINCOMP_K{1,i} = [Tincon,Pincon(:,6)];
 
end 
figure('Name','Non Competitiva analisi di SENSITIVITA')

hold on
for i = 1:length(k3new)
 

 subplot(2,4,1)                                   % S NO COMP
 plot(sNO_COMP_K{1,i}(:,1),sNO_COMP_K{1,i}(:,2))
 title('substrato [S]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,2)                                   % E NO COMP
 plot(eNO_COMP_K{1,i}(:,1),eNO_COMP_K{1,i}(:,2))
 title('enzima [E]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,3)                                   % ES NO COMP
 plot(cNO_COMP_K{1,i}(:,1),cNO_COMP_K{1,i}(:,2))
 title('complesso [ES]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,4)                                   % P NO COMP
 plot(pNO_COMP_K{1,i}(:,1),pNO_COMP_K{1,i}(:,2))
 title('prodotto [P]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,5)                                   % I NO COMP
 plot(iNO_COMP_K{1,i}(:,1),iNO_COMP_K{1,i}(:,2))
 title('inibitore [I]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,6)                                   % EI NO COMP
 plot(eiNO_COMP_K{1,i}(:,1),eiNO_COMP_K{1,i}(:,2))
 title('complesso 2 [EI]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on
 
 subplot(2,4,7)                                   % EIS NO COMP
 plot(eisNO_COMP_K{1,i}(:,1),eisNO_COMP_K{1,i}(:,2))
 title('prodotto 2 [EIS]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on
 

 legend({'k3=25*10^2','k3=50*10^2','k3=75*10^2','k3=150*10^2','k3=250*10^2','k3=500*10^2'}, ...
 'Position',[0.747,0.109,0.156,0.344]) 
 

end 


figure('Name','Competitiva analisi di SENSITIVITA')

hold on
for i = 1:length(k3new)
 

 subplot(2,4,1)                                   % S COMP
 plot(sCOMP_K{1,i}(:,1),sCOMP_K{1,i}(:,2))
 title('substrato [S]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,2)                                   % E COMP
 plot(eCOMP_K{1,i}(:,1),eCOMP_K{1,i}(:,2))
 title('enzima [E]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,3)                                   % ES COMP
 plot(cCOMP_K{1,i}(:,1),cCOMP_K{1,i}(:,2))
 title('complesso [ES]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,4)                                   % P COMP
 plot(pCOMP_K{1,i}(:,1),pCOMP_K{1,i}(:,2))
 title('prodotto [P]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,5)                                   % I COMP
 plot(iCOMP_K{1,i}(:,1),iCOMP_K{1,i}(:,2))
 title('inibitore [I]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,6)                                   % EI COMP
 plot(eiCOMP_K{1,i}(:,1),eiCOMP_K{1,i}(:,2))
 title('complesso 2 [EI]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on
 

 legend({'k3=25*10^2','k3=50*10^2','k3=75*10^2','k3=150*10^2','k3=250*10^2','k3=500*10^2'}, ...
 'Position',[0.747,0.109,0.156,0.344]) 

end 
figure('Name','Incompetitiva analisi di SENSITIVITA')
% fig 22
% ATTENZIONE:
% non è che non mi stampa le righe è che sono sovrapposte
% la variazione di k1 non influenza l'andamento
hold on
for i = 1:length(k3new)
 

 subplot(2,4,1)                                   % S INCOMP
 plot(sINCOMP_K{1,i}(:,1),sINCOMP_K{1,i}(:,2))
 title('substrato [S]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,2)                                   % E INCOMP
 plot(eINCOMP_K{1,i}(:,1),eINCOMP_K{1,i}(:,2))
 title('enzima [E]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,3)                                   % ES INCOMP
 plot(cINCOMP_K{1,i}(:,1),cINCOMP_K{1,i}(:,2))
 title('complesso [ES]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,4)                                   % P INCOMP
 plot(pINCOMP_K{1,i}(:,1),pINCOMP_K{1,i}(:,2))
 title('prodotto [P]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,5)                                   % I INCOMP
 plot(iINCOMP_K{1,i}(:,1),iINCOMP_K{1,i}(:,2))
 title('inibitore [I]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 hold on

 subplot(2,4,6)                                   % EIS INCOMP
 plot(eisINCOMP_K{1,i}(:,1),eisINCOMP_K{1,i}(:,2))
 title('prodotto 2 [EIS]')
 xlabel('Tempo [sec]')
 ylabel('Concentrazione [mM]')
 

 hold on
 legend({'k3=25*10^2','k3=50*10^2','k3=75*10^2','k3=150*10^2','k3=250*10^2','k3=500*10^2'}, ...
 'Position',[0.747,0.109,0.156,0.344]) 

 
end 
%% V0 SU S0

SSS=[0:1:1000];
KM =(k_1e+k2e)/k1e;
ki = k_3in/k3in;
ks = (k_1in+k2in)/k1in;
VMAX = k2e.*e0.*SSS./(SSS+KM);

figure('Name','confronto')

subplot(4,2,1)                                          % STANDARD
plot(SSS,VMAX)

title('normale')
xlabel('S0')
ylabel('V0')

subplot(4,2,2)
plot(1./SSS,1./VMAX)
title('normale')
xlabel('1/S0')
ylabel('1/V0')

I0 = [0.1 1 10];

for i=1:3
 VMAX_COMP = (SSS./ks)./(1+(SSS./ks)+(I0(i)./ki));
 VMAX_INC = (SSS./ks)./(1+(SSS./ks)+((SSS./ks)*(I0(i)./ki)));
 VMAX_NC = (SSS./ks)./(1+(SSS./ks)+(I0(i)./ki)+((SSS./ks)*(I0(i)./ki)));
 
 

 subplot(4,2,3)                              % NO COMP
 plot(SSS,VMAX_NC)
 hold on
 legend('I=0.1','I=1','I=10')
 title('non competitiva')
 xlabel('S0')
 ylabel('V0')


 subplot(4,2,4)                              % NO COMP Lineweaver-Burk Plot
 plot(1./SSS,1./VMAX_NC)
 hold on
 legend('I=0.1','I=1','I=10')
 title('non competitiva')
 xlabel('1/S0')
 ylabel('1/V0')



 subplot(4,2,5)                              % INCOMP
 plot(SSS,VMAX_INC)
 hold on
 legend('I=0.1','I=1','I=10')
 title('incompetitiva')
 xlabel('S0')
 ylabel('V0')


 subplot(4,2,6)                              % INCOMP Lineweaver-Burk Plot
 plot(1./SSS,1./VMAX_INC)
 hold on
 legend('I=0.1','I=1','I=10')
 title('incompetitiva')
 xlabel('1/S0')
 ylabel('1/V0')



 subplot(4,2,7)                              % COMP
 plot(SSS,VMAX_COMP)
 hold on
 legend('I=0.1','I=1','I=10')
 title('competitiva')
 xlabel('S0')
 ylabel('V0')


 subplot(4,2,8)                              % COMP Lineweaver-Burk Plot
 plot(1./SSS,1./VMAX_COMP)
 hold on
 legend('I=0.1','I=1','I=10')
 title('competitiva')
 xlabel('1/S0')
 ylabel('1/V0')
end 