%% ES 3

%%
clear 
close all
V1 = 3.29;     %[l]
k01 = 6.54e-2; %[1/min]
k12 = 5.68e-2; %[1/min]
k21 = 7.16e-2; %[1/min]
D_es = 49650;  %[pmol]


A1 = [-(k01+k21), k12; k21, -k12];
B1 = [D_es;0];
C1 = [1/V1, 0];
D = 0;

Sistema1 = ss(A1,B1,C1,D);

[Concentrazione, Tempo, Q] = impulse(Sistema1);


figure(1)
subplot(3,1,1)
plot(Tempo,Concentrazione)

title('CONCENTRAZIONE compartimento 1');
xlabel('t [min]');
ylabel('Concentrazione [pmol/l]');

% GRAFICO log(Concentrazione)

subplot(3,1,2)

% CONVERTO SU SCALA LOGARITMICA

semilogy(Tempo,Concentrazione)

title('LOGARITMO DELLA CONCENTRAZIONE');
xlabel('t [min]');
ylabel('Concentrazione [pmol/l]');

% DUE PENDENZE ---> DUE COSTANTI

%GRAFICO Q1

subplot(3,1,3)

plot(Tempo,Q)
title('QUANTITA compartimento 1 E 2');
xlabel('tempo [min]');
ylabel('Quantità [pmol]');

%% A)

transfer_function = tf(Sistema1);

autovalori = eig(A1);

tau1 = -1/autovalori(1);  %tau veloce
tau2 = -1/autovalori(2);  %tau lenta

tempo_eliminazione=5*tau2;  %[min]


% IL TEMPO DI ELIMINAZIONE DIPENDERA' SEMPRE DALLA COSTANTE DI ELIMINAZIONE
% PIU' LENTA


%% B)

% LA CONCENTRAZIONE MASSIMA MISURABILE DIPENDE DA:
% INGRESSO ESOGENO
% COSTANTI DI SCAMBIO E ELIMINAZIONE
% VOLUME DEI COMPARTIMENTI

Concentrazione_MAX = max(Concentrazione);

%% C)

% COME LA CONCENTRAZIONE MASSIMA, L'AREA SOTTO LA CURVA DIPENDERA' SEMPRE
% DA TUTTI I PARAMETRI ELENCATI IN PRECEDENZA

AUC = trapz(Tempo,Concentrazione);

AUC2 = D_es/(V1*k01);

%% D)

CL = V1*k01;

% DIPENDERA' SOLO DAL VOLUME 1 E DALLA COSTANTE DI ELIMINAZIONE k01,
% PERCHE' E' L'UNICO COMPARTIMENTO IN CUI VIENE ELIMINATO



%% E)

k21_E = [k21/100;...
         k21/80;...
         k21/60;...
         k21/40;...
         k21/20;...
         k21*1;...
         k21*10;...
         k21*20;...
         k21*50;...
         k21*100];

figure(3)

for i=1:10
    sys_cycle=ss([-(k01+k21_E(i)) k12; k21_E(i) -k12],[D_es; 0],[1/V1 0],0);
    [Y_cycle,T_cycle]=impulse(sys_cycle);
    subplot(2,1,1)
    plot(T_cycle,Y_cycle), grid on
    title('CONCENTRAZIONE cinetica primo compartimento');
    legend('k21/100','k21/80','k21/60','k21/40','k21/20','k21*1','k21*10','k21*20','k21*50','k21*100');
    hold on
    xlabel('t [min]');
    ylabel('Concentrazione [pmol/l]');
    
    subplot(2,1,2)
    plot(T_cycle,Y_cycle.*V1), grid on
    title('QUANTITA di farmaco primo compartimento');
    hold on
    xlabel('t [min]');
    ylabel('Q [pmol]');
end

% PIU' K21 TENDERA' A 0 IL CONTRIBUTO RIMANE ESCLUSIVAMENTE QUELLO DELLA
% k01 PRESENTE NEL PRIMO COMPARTIMENTO


%% DATI REALISTICI 2022
% set realistico di dati con 7 punti della curva
t_camp2 = linspace(Tempo(1),Tempo(end),7);
indici2 = floor(linspace(1,346,7));
c_camp2 = Concentrazione(indici2);

AUC_camp2 = trapz(t_camp2,c_camp2);


figure(4)

plot(Tempo,Concentrazione)
hold on

xlabel('tempo [ore]')
ylabel('concentrazione [mg/l]')
hold on 
scatter(t_camp2,c_camp2);
