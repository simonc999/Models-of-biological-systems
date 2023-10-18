%% ES 1
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


%% A)B) TEMPO DI ELIMINAZIONE E DI EMIVITA

tau1 = 1/k01; % [h]

% IL TEMPO DI ELIMINAZIONE SARA' 5*TAU

tempo_eliminazione = 5*tau1;

% IL TEMPO DI ELIMINAZIONE DIPENDE DI FATTO DALLA COSTANTE DI ELIMINAZIONE
% k01, PIU' SARA' ELEVATA PIU' IL TEMPO DI ELIMINAZIONE SARA' RAPIDO

% B) TEMPO DI EMIVITA

% TEMPO PER CUI LA CONCENTRAZIONE E' LA META' DELLA CONCENTRAZIONE INIZIALE

tempo_emivita = log(2)*tau1;

% AVREMO UN VALORE PARI AL TEMPO DI EMIVITA NELLA POSIZIONE 16 DEL VETTORE
% DEL TEMPO

% POSSIAMO INFATTI NOTARE CHE NELLA POSIZIONE 16 DELLE QUANTITA' E DELLE
% CONCENTRAZIONI I DUE PARAMETRI SONO LA META DEI VALORI INZIZIALI

C_EMIVITA = Concentrazione(16);   % META' DI 100 mg/l

Q_EMIVITA = Q1(16);               % META' DI 500 mg


%% C) CONCENTRAZIONE MASSIMA MISURABILE

% DIPENDE DALL'INGRESSO ESOGENO

% MASSIMA ALL'ISTANTE 0 CONSIDERANDO LA FUNZIONE ESPONENZIALE

% MASSIMA ALL'ISTANTE 0 CONSIDERANDO UNA SOMMINISTRAZIONE ISTANTANEA PARI A
% D_es

Concentrazione_MAX1 = Concentrazione(1);    %[mg/l]
Concentrazione_MAX2 = max(Concentrazione);  %[mg/l]

%% D)
% L'AREA SOTTO LA CURVA DIPENDE DI FATTO DA TUTTI I PARAMETRI PRESENTI NEL
% SISTEMA

% QUINDI DA D_es, V1 e k01

AUC_1 = trapz(Tempo,Concentrazione);

AUC_2 = D_es/(k01*V1);   % SEGUENDO LA DIMOSTRAZIONE DELL'ESPONENZIALE (FOGLIO)

AUC_3 = Concentrazione_MAX1/k01;

%% E) 

AUMC=trapz(Tempo,Tempo.*Concentrazione);

tempo_medio_residenza=V1/(V1*k01);
tempo_medio_residenza2=AUMC/AUC_1;  % MRT

%% F)

CL_tot=k01*V1;

%% G)
%% ------------------------------------------------------------------------
%% ------------------------------TEST--------------------------------------
%% -----------------------------VOLUME-------------------------------------
%% ------------------------------------------------------------------------

% TEST VOLUME--QUADRUPLICATO E DIVISO PER 4

% V1=4*V1

V1_GV1 = 20;    %[l]

% V1=0.25*V1

V1_GV2 = 1.25;  %[l]

% IMPOSTO IL SISTEMA PER GV1

A1_GV1 = -k01;
B1_GV1 = D_es;
C1_GV1 = 1/V1_GV1;
D1_GV1 = 0;

Sistema1_GV1 = ss(A1_GV1,B1_GV1,C1_GV1,D1_GV1);

[Concentrazione_GV1,Tempo_GV1] = impulse(Sistema1_GV1);

Q1_GV1 = Concentrazione_GV1*V1_GV1;

% IMPOSTO IL SISTEMA PER GV2

A1_GV2 = -k01;
B1_GV2 = D_es;
C1_GV2 = 1/V1_GV2;
D1_GV2 = 0;

Sistema1_GV2 = ss(A1_GV2,B1_GV2,C1_GV2,D1_GV2);

[Concentrazione_GV2,Tempo_GV2] = impulse(Sistema1_GV2);

Q1_GV2 = Concentrazione_GV2*V1_GV2;

figure(3)

%-----------------------------GV1------------------------------------

%  GRAFICO Concentrazione GV1 

subplot(3,2,1)
plot(Tempo_GV1,Concentrazione_GV1)
title('Concentrazione GV1 4V');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO LOG(Concentrazione_GV1)

subplot(3,2,3)
semilogy(Tempo_GV1,Concentrazione_GV1)
title('Concentrazione GV1 4V su scala logaritmica');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO Q_GV1

subplot(3,2,5)
plot(Tempo_GV1,Q1_GV1,"Color",[1 0 0])
title('QUANTITA GV1 4V');
xlabel('Tempo [h]');
ylabel('Quantità [mg]');

%-----------------------------GV2------------------------------------

%  GRAFICO Concentrazione GV2 

subplot(3,2,2)
plot(Tempo_GV2,Concentrazione_GV2)
title('Concentrazione GV2 0.25V');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO LOG(Concentrazione_GV2)

subplot(3,2,4)
semilogy(Tempo_GV2,Concentrazione_GV2)
title('Concentrazione GV2 0.25V su scala logaritmica');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO Q_GV2

subplot(3,2,6)
plot(Tempo_GV2,Q1_GV2,"Color",[1 0 0])
title('QUANTITA GV2 0.25V');
xlabel('Tempo [h]');
ylabel('Quantità [mg]');


% OSSERVO CHE ALL'AUMENTARE DEL VOLUME LA CONCENTRAZIONE DECRESCE MENTRE
% L'ANDAMENTO DELLA QUANTITA RIMANE COSTANTE. LE DUE GRANDEZZE DATE DALLA
% SONO INVERSAMENTE PROPORZIONALI, E SONO LEGATE DALLA RELAZIONE DI
% QUANTITA' DI FARMACO CALCOLATA CHE DI FATTO RIMANE COSTANTE IN TUTTI I
% CASI

% DI FATTO IL FLUSSO Q PUNTO INTERESSA LA QUANTITA' NEL TEMPO, QUINDI
% DIPENDERA' SOLTANTO DA D E DA k01

%% ------------------------------------------------------------------------
%% ------------------------------TEST--------------------------------------
%% ---------------------COSTANTE DI ELIMINAZIONE---------------------------
%% ------------------------------------------------------------------------
% IN MODO ANALOGO RADDOPPIO E DIMEZZO LA COSTANTE DI ELIMINAZIONE

k01_GK1=2.4; %[h^-1]

k01_GK2=0.6; %[h^-1]

% MODIFICANDO LA COSTANTE DI ELIMINAZIONE SI PUO' NOTARE CHE, ESSENDO
% L'UNITA' DI MISURA DI k01 ore^-1, PIU' SI INCREMENTA TALE VALORE PIU' IL
% FARMACO VERRA' ELIMINATO VELOCEMENTE, IN MODO CONTRARIO, PIU' VERRA
% DIMINUITA E PIU' IL FARMACO VERRA' ELIMINATO LENTAMENTE

% NOTIAMO INFATTI CHE TAU E' INVERSAMENTE PROPORZIONALE A k01

% IL TEMPO DI ELIMINAZIONE PER IL SISTEMA LA CUI COSTANTE E' k01_GK1 SARA'
% LA META' DI QUELLO ORIGINALE, MENTRE PER IL SISTEMA LA CUI COSTANTE E'
% k01_GK2 SARA' RADDOPPIATO


TAU_GK1=1/k01_GK1;

TEMPO_ELIMINAZIONE_GK1=5*TAU_GK1;

TAU_GK2=1/k01_GK2;

TEMPO_ELIMINAZIONE_GK2=5*TAU_GK2;


% IMPOSTO IL SISTEMA PER CONCENTRAZIONE RADDOPPIATA

A1_GK1 = -k01_GK1;
B1_GK1 = D_es;
C1_GK1 = 1/V1;
D1_GK1 = 0;

Sistema1_GK1 = ss(A1_GK1,B1_GK1,C1_GK1,D1_GK1);

[Concentrazione_GK1,Tempo_GK1] = impulse(Sistema1_GK1);

Q1_GK1 = Concentrazione_GK1*V1;


% IMPOSTO IL SISTEMA PER CONCENTRAZIONE DIMEZZATA

A1_GK2 = -k01_GK2;
B1_GK2 = D_es;
C1_GK2 = 1/V1;
D1_GK2 = 0;

Sistema1_GK2 = ss(A1_GK2,B1_GK2,C1_GK2,D1_GK2);

[Concentrazione_GK2,Tempo_GK2] = impulse(Sistema1_GK2);

Q1_GK2 = Concentrazione_GK2*V1;


figure(4)

%-----------------------------GK1------------------------------------

%  GRAFICO Concentrazione GK1 

subplot(3,2,1)
plot(Tempo_GK1,Concentrazione_GK1)
title('Concentrazione GK1 2k01');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');

xline(TEMPO_ELIMINAZIONE_GK1,'--r',LineWidth=2);
% GRAFICO LOG(Concentrazione_GK1)

subplot(3,2,3)
semilogy(Tempo_GK1,Concentrazione_GK1)
title('Concentrazione GK1 2k01 su scala logaritmica');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO Q_GK1

subplot(3,2,5)
plot(Tempo_GK1,Q1_GK1,"Color",[1 0 0])
title('QUANTITA GK1 2k01');
xlabel('Tempo [h]');
ylabel('Quantità [mg]');

%-----------------------------GK2------------------------------------

%  GRAFICO Concentrazione GK2 

subplot(3,2,2)
plot(Tempo_GK2,Concentrazione_GK2)
title('Concentrazione GK2 0.5k01');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');
xline(TEMPO_ELIMINAZIONE_GK2,'--r',LineWidth=2);
% GRAFICO LOG(Concentrazione_GK2)

subplot(3,2,4)
semilogy(Tempo_GK2,Concentrazione_GK2)
title('Concentrazione Gk2 0.5k01 su scala logaritmica');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO Q_GK2

subplot(3,2,6)
plot(Tempo_GK2,Q1_GK2,"Color",[1 0 0])
title('QUANTITA GK2 0.5k01');
xlabel('Tempo [h]');
ylabel('Quantità [mg]');




%% ------------------------------------------------------------------------
%% ------------------------------TEST--------------------------------------
%% ------------------------INGRESSO ESOGENO--------------------------------
%% ------------------------------------------------------------------------

D_es_GD1=100;  %[mg]

D_es_GD2=1000; %[mg]

% IMPOSTO IL SISTEMA GD1 CON INGRESSO 1/5 QUELLO ORIGINARIO

A1_GD1 = -k01;
B1_GD1 = D_es_GD1;
C1_GD1 = 1/V1;
D1_GD1 = 0;

Sistema1_GD1 = ss(A1_GD1,B1_GD1,C1_GD1,D1_GD1);

[Concentrazione_GD1,Tempo_GD1] = impulse(Sistema1_GD1);

Q1_GD1 = Concentrazione_GD1*V1;

% IMPOSTO IL SISTEMA GD2 CON INGRESSO 2 VOLTE QUELLO ORIGINARIO

A1_GD2 = -k01;
B1_GD2 = D_es_GD2;
C1_GD2 = 1/V1;
D1_GD2 = 0;

Sistema1_GD2 = ss(A1_GD2,B1_GD2,C1_GD2,D1_GD2);

[Concentrazione_GD2,Tempo_GD2] = impulse(Sistema1_GD2);

Q1_GD2 = Concentrazione_GD2*V1;


figure(5)
%-----------------------------GD1------------------------------------

%  GRAFICO Concentrazione GD1 

subplot(3,2,1)
plot(Tempo_GD1,Concentrazione_GD1)
title('Concentrazione GD1 0.2D_es');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO LOG(Concentrazione_GK1)

subplot(3,2,3)
semilogy(Tempo_GD1,Concentrazione_GD1)
title('Concentrazione GD1 0.2D_es su scala logaritmica');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO Q_GK1

subplot(3,2,5)
plot(Tempo_GD1,Q1_GD1,"Color",[1 0 0])
title('QUANTITA GD1 0.2D_es');
xlabel('Tempo [h]');
ylabel('Quantità [mg]');

%-----------------------------GD2------------------------------------

%  GRAFICO Concentrazione GD2 

subplot(3,2,2)
plot(Tempo_GD2,Concentrazione_GD2)
title('Concentrazione GD2 2D_es');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO LOG(Concentrazione_GD2)

subplot(3,2,4)
semilogy(Tempo_GD2,Concentrazione_GD2)
title('Concentrazione GD2 2D_es su scala logaritmica');
xlabel('Tempo [h]');
ylabel('Concentrazione [mg/l]');

% GRAFICO Q_GD2

subplot(3,2,6)
plot(Tempo_GD2,Q1_GD2,"Color",[1 0 0])
title('QUANTITA GD2 2D_es');
xlabel('Tempo [h]');
ylabel('Quantità [mg]');


% A VOLUME COSTANTE E PARAMETRO K DI ELIMINAZIONE COSTANTE AVRO' UNA
% CONCENTRAZIONE MINORE O MAGGIORE PROPORZIONALE DIRETTAMENTE ALLA
% QUANTITA' DI FARMACO IN INGRESSO ESOGENO


% SE K RIMANE COSTANTE ANCHE LA K DI ELIMINAZIONE RIMANE COSTANTE 


%% DATI REALISTICI 2022

% uso un set di dati realistici
% considero al massimo 7 prelievi 
% attuo un campionamento

t_camp = linspace(Tempo(1),Tempo(end),7);
index = linspace(1,127,7);
c_camp = Concentrazione(index);

% noto che se uno questo set sbaglio ma di poco rispetto alla AUC calcolata
% sull'intera curva non campionata

AUC_camp = trapz(t_camp,c_camp);
CL_realistica = D_es/AUC_camp;

figure(6)

plot(Tempo,Concentrazione)
xlabel('tempo [ore]')
ylabel('concentrazione [mg/l]')
hold on 
scatter(t_camp,c_camp);
