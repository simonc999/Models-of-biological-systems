%% ES 2

%% MODELLO A UN COMPARTIMENTO CON ASSORBIMENTO
clear
close all
V2 = 5;         %[l]

k01 = 1.2;      %[h^-1]
k02 = 1.2;      %[h^-1]
k21 = 2.2;      %[h^-1]

D_es = 500;     %[mg] bolo per via orale 

A2 = [-(k01+k21),  0; k21,  -k02];
B2 = [D_es;0];
C2 = [0, 1/V2];
D  = 0;

Sistema2 = ss(A2,B2,C2,D);

[Concentrazione_2, Tempo, Quantita_12] = impulse(Sistema2);


%  GRAFICO Concentrazione

figure(1)
subplot(3,1,1)
plot(Tempo,Concentrazione_2)

title('CONCENTRAZIONE compartimento 2 (cinetica mista)');
xlabel('t [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO log(Concentrazione)

subplot(3,1,2)

% CONVERTO SU SCALA LOGARITMICA

semilogy(Tempo,Concentrazione_2)

title('LOGARITMO DELLA CONCENTRAZIONE');
xlabel('t [h]');
ylabel('Concentrazione [mg/l]');

% DUE PENDENZE ---> DUE COSTANTI DI ELIMINAZIONE

%GRAFICO Q1

subplot(3,1,3)


plot(Tempo,Quantita_12)
title('QUANTITA compartimento 1 e 2');
xlabel('tempo [h]');
ylabel('Quantità [mg]');


%% A)

autovalori=eig(A2);

tau1 = -1/autovalori(1); % OPPURE 1/1.2       (k02)     [h] 
tau2 = -1/autovalori(2); % OPPURE 1/(1.2+2.2) (k21+k01) [h]

tempo_eliminazione = 5*tau1;
tempo_assorbimento = 5*tau2;

% CONSIDERO SEMPRE LA PIU' LENTA

% F(s)

transfer_function = tf(Sistema2);

[R, P, K] = residue(transfer_function.Numerator{1}, transfer_function.Denominator{1});

%% B)

CONCENTRAZIONE_2_INIT = Concentrazione_2(1);

CONCENTRAZIONE_2_MASSIMA = max(Concentrazione_2);

% LA CONCENTRAZIONE MASSIMA MISURABILE DIPENDERA' DAI VOLUMI DEI DUE
% COMPARTIMENTI, DAL BOLO SOMMINISTRATO E DALLE COSTANTI DI ASSORBIMENTO E
% DI ELIMINAZIONE

%% C)

AUC = trapz(Tempo,Concentrazione_2);


%% D)

% DALLA TEORIA

% La frazione biodisponibile (o biodisponibilità) è la frazione della
% dose somministrata che raggiunge effettivamente il sistema
% circolatorio. Per definizione quando un farmaco è somminsitrato
% per via endovenosa è biodisponibile al 100%. Quando il farmaco è
% somministrato per altra via la biodisponibilià generalemente
% decresce a causa dell’incompleto assorbimento

F = k21/(k01+k21);

% LA DOSE BIODISPONIBILE SI OTTERRA' MOLTIPLICANDO LA FRAZIONE
% BIODISPONIBILE PER LA DOSE SOMMINISTRATA

DOSE_F = F*D_es;

%% E)

k01_E = [k01/100;...
         k01/80;...
         k01/60;...
         k01/40;...
         k01/20;...
         k01*1;...
         k01*10;...
         k01*20;...
         k01*50;...
         k01*100];

figure(3)

for i=1:numel(k01_E)

    A2_cy = [-(k01_E(i)+k21),  0; k21,  -k02];
    B2_cy = [D_es;0];
    C2_cy = [0, 1/V2];
    D_cy  = 0;
    
    sys_cycle = ss(A2_cy,B2_cy,C2_cy,D_cy);

    [Y_cycle,T_cycle] = impulse(sys_cycle);
    subplot(2,1,1)
    plot(T_cycle,Y_cycle), grid on
    title('CONCENTRAZIONE cinetica secondo compartimento');
    legend('k01/100','k01/80','k01/60','k01/40','k01/20','k01*1','k01*10','k01*20','k01*50','k01*100');
    hold on
    xlabel('t [h]');
    ylabel('Concentrazione [mg/l]');
    
    subplot(2,1,2)
    plot(T_cycle,Y_cycle*V2), grid on
    title('QUANTITA di farmaco primo e secondo compartimento');
    hold on
    xlabel('t [h]');
    ylabel('Q [mg]');
end
%% F)

% PER CALCOLARE LA FRAZIONE DI BIODISPONIBILITA' CON UN APPROCCIO NON 
% COMPARTIMENTALE ABBIAMO BISOGNO DEI PARAMETRI AUC E D DELL'ESERCIZIO 1 , 
% CHE SUPPONEVA L'ANALISI COMPARTIMENTALE PER UN INGRESSO DI TIPO 
% ENDOVENOSO, QUINDI SENZA RIASSORBIMENTO

AUC_ev = 83.333330;

AUMC_ev = 68.0054;

AUMC = trapz(Tempo,Concentrazione_2.*Tempo);

% LA DOSE ORALE SARA' PARI ALLA DOSE PER LA FRAZIONE BIODISPONIBILE

F_non_compartimentale = AUC*D_es/(AUC_ev*D_es);

%% G)

% 1/ka = MAT    MAT = MRTorale - MRTev

MRTev = AUMC_ev/AUC_ev;

MRTorale = AUMC/AUC;

MAT = MRTorale - MRTev;

ka = 1/MAT;   % COSTANTE DI ASSORBIMENTO O VELOCITA DI ASSORBIMENTO


%% DATI REALISTICI 2022

t_camp1 = linspace(Tempo(1),Tempo(end),7);
indici1 = floor(linspace(1,230,7));
c_camp1 = Concentrazione_2(indici1);

AUC_camp1 = trapz(t_camp1,c_camp1);
CL_realistica1 = D_es/AUC_camp1;

figure(7)
plot(Tempo,Concentrazione_2)
xlabel('tempo [ore]')
ylabel('concentrazione [mg/l]')
hold on 
scatter(t_camp1,c_camp1);