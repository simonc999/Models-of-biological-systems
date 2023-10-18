%% ES 5
clear 
close all
% MODELLO MONOCOMPARTIMENTALE

% SOMMINISTRAZIONE BOLO

% INGRESSO ESOGENO

D_es = 500;        %[mg]

% VOLUME COMPARTIMENTO 1

V1 = 5;          %[l]

%  COSTANTE DI ELIMINAZIONE k01

k01 = 1.2;      %[h^-1]

% IMPOSTO IL SISTEMA

A1 = -k01;
B1 = D_es;
C1 = 1/V1;
D1 = 0;

Sistema1 =ss (A1,B1,C1,D1);

[Concentrazione,Tempo] = impulse(Sistema1);  % [mg/l] , [h]

Q1=V1.*Concentrazione;

%  GRAFICO Concentrazione

figure(1)
subplot(3,1,1)
plot(Tempo,Concentrazione)

title('CONCENTRAZIONE compartimento 1 (cinetica primo ordine monocompartimentale)');
xlabel('t [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO log(Concentrazione)

subplot(3,1,2)

% CONVERTO SU SCALA LOGARITMICA

semilogy(Tempo,Concentrazione)

title('LOGARITMO DELLA CONCENTRAZIONE');
xlabel('t [h]');
ylabel('Concentrazione [mg/l]');

% UNA SOLA PENDENZA ---> UNA SOLA COSTANTE (k01) DI ELIMINAZIONE

%GRAFICO Q1

subplot(3,1,3)


plot(Tempo,Q1,"Color",[1 0 0])
title('QUANTITA compartimento 1 (cinetica primo ordine monocompartimentale)');
xlabel('tempo [h]');
ylabel('Quantità [mg]');

%% 
ampiezzaImpulsi=(max(Concentrazione)-min(Concentrazione))/(C1*B1); %quantita di farmaco da iniettaare, le ottengo con formule di lagrange
periodo=6;
numeroImpulsi=12;
durata=3*24;

t = linspace(0,periodo,durata);
x = 0;
Y = [];
X = [];
T = [];
for i = 1:numeroImpulsi
   [yl,tl,xl] = initial(Sistema1,x,t);
   [yf,tf,xf] = impulse(Sistema1*ampiezzaImpulsi,t);
   Y = [Y; yl+yf];
   X = [X; xl+xf];
   T = [T t+periodo*(i-1)];
   x = xf(end,:) + xl(end,:);
end

figure(3)
subplot(2,1,1)
plot(T,Y), grid on
title('ANDAMENTO CONCENTRAZIONE TRENO DI IMPULSI')
xlabel('tempo [h]');
ylabel('Concentrazione [mg/l]');

subplot(2,1,2)
plot(T,X), grid on
title('RISPOSTA DELLA QUANTITA')
xlabel('tempo [h]');
ylabel('Quantità [mg]');


%% A)

tau1 = 1/k01; % [h]

% IL TEMPO DI ELIMINAZIONE SARA' 5*TAU

tempo_eliminazione = 5*tau1;

% IL TEMPO DI ELIMINAZIONE DIPENDE DI FATTO DALLA COSTANTE DI ELIMINAZIONE
% k01, PIU' SARA' ELEVATA PIU' IL TEMPO DI ELIMINAZIONE SARA' RAPIDO

%% B)

CSS = D_es/(V1*k01*periodo);

%% C) FOGLIO

%% D)


% Somministrazione ORALE

V2 = 5;         %[l]

k01 = 1.2;      %[h^-1]
k02 = 1.2;      %[h^-1]
k21 = 2.2;      %[h^-1]

D_es = 500;     %[mg] bolo per via orale 

A2 = [-(k01+k21),  0; k21,  -k02];
B2 = [1;0];
C2 = [0, 1/V2];
D  = 0;

% usando il modello compartimentale del es 2
% varia la matrice B quindi rifaccio il sistema
% A1 C1 D1 rimagono le stesse e quindi le riciclo
% la durata totale è invariata = 7 gg
% la concentrazione rimane 50 mg/l
% con una variazione 10%


FRAZIONE_K= (k21/(k02-k01-k21));

% uno dei due non è esponenziale quindi non posso fare nulla 
% vafo a tentativi

SYSor = ss(A2,B2,C2,D);
durata_tot = 24*7;        % 7 giorni
AI_orale = 50/0.647;          % a tentativi
Tdecad_orale = 0.17;        % a tentativi
f_orale = 1/Tdecad_orale;
NI_orale = durata_tot/Tdecad_orale;

tor = linspace(0,Tdecad_orale,durata_tot);
xor = [0 0];
Yor = [];
Xor = [];
Tor = [];
for i = 1:NI_orale
   [ylor,tlor,xlor] = initial(SYSor,xor,tor);
   [yfor,tfor,xfor] = impulse(SYSor*AI_orale,tor);
   Yor = [Yor; ylor+yfor];
   Xor = [Xor; xlor+xfor];
   Tor = [Tor tor+Tdecad_orale*(i-1)];
   xor = xfor(end,:) + xlor(end,:);
end

figure(15)
plot(Tor,Yor)
grid on
title('Risposta al treno di impulsi')









