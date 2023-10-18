%% ESERCITAZIONE 3 MODELLI

% ATTENZIONE a fine script viene stmpato un REPORT riassuntivo
%            titolo del report = 'eserc3modelli.txt'
clc
clear 
close all

%%% risposta ad un bolo di C-peptide sintetico in 7 pz.
% CINETICA LINEARE

% C PEP PRODUCE PRO INSULINA  POI PRE FINO AL PANCREAS DOVE VIENE PRODOTTA
% L'INSULINA

% C PEP HA CINETICA LINEARE A DIFFERENZA DELL'INSULINA 
% L'INSULINA SECRETA NEL PANCREAS E ASSORBIMENTO NEL FEGATO
% C PEP E' COMPONENTE INDICATORE DI SECREZIONE DELL'INSULINA
% SECRETI NELLA STESSA QUANTITA' 
% NEL PLASMA CONC DI C MAGGIORE
% MISURIAMO IL C PER DETERMINARE LA PRODUZIONE DELL'INSULINA DA PARTE DELLE
% CELLULE BETA DEL PANCREAS

% MISURO IL C PERCHE' DA 100 DI C SICURO MENO DI 100 DI INSULINA (PROCESSI
% INTERMEDI)

% C SIGNIFICATO COMPLEMNETARE (DELL'INSULINA) POI INSULINA VA VIA

% IN SEGUITO ALLA SOMMINISTRAZIONE DELLA SOMATOSTRATINA ( INIBITORE DEL C
% ENDOGENO) NOI IMMETTIAMO UNA DOSE DI C SINTETICO IN MODO DA OSSERVARNE LA
% CINETICA


% MINUTO 0 MISURO C ENDOGENO CHE E' RIMASTO A SEGUITO DI INIBIZIONE


S1 = table2array(readtable('DatiCPsog1.dat'));
S2 = table2array(readtable('DatiCPsog2.dat'));
S3 = table2array(readtable('DatiCPsog3.dat'));
S4 = table2array(readtable('DatiCPsog4.dat'));
S5 = table2array(readtable('DatiCPsog5.dat'));
S6 = table2array(readtable('DatiCPsog6.dat'));
S7 = table2array(readtable('DatiCPsog7.dat'));

figure(1)
subplot(2,1,1)
plot(S1(:,1),S1(:,2))
hold on
plot(S2(:,1),S2(:,2))
hold on
plot(S3(:,1),S3(:,2))
hold on
plot(S4(:,1),S4(:,2),'LineWidth',2)
hold on
plot(S5(:,1),S5(:,2))
hold on
plot(S6(:,1),S6(:,2))
hold on
plot(S7(:,1),S7(:,2))
hold on
legend('1','2','3','4','5','6','7')
xlabel('Tempo [minuti]')
ylabel('Concentrazione [pmol/ml]')

subplot(2,1,2)
semilogy(S1(:,2))
hold on
semilogy(S2(:,2))
hold on
semilogy(S3(:,2))
hold on
semilogy(S4(:,2),'LineWidth',2)
hold on
semilogy(S5(:,2))
hold on
semilogy(S6(:,2))
hold on
semilogy(S7(:,2))
hold on
legend('1','2','3','4','5','6','7')

% CONSIDERO QUARTO SOGGETTO (LA CURVA STA IN MEZZO NELLA SCALA SEMILOG)

TEMPO = S4(:,1);                       % [minuti]
CONCENTRAZIONE = S4(:,2);              % [pmol/ml]
BASALE = CONCENTRAZIONE(TEMPO==0);


% SOTTRAGGO la concentrazione BASALE
%   Ho un sistema lineare: posso sfruttare la sovrapposizone degli effetti.
%   Voglio esaminare la risposta al bolo quindi elimino ciò che era già 
%   presente nell'organismo, studio il sistema alle differenze

% ELIMINO la misura effettuata al MINUTO 1
%   essa deve essere trascurata nell identificazione del modello in quanto 
%   ci sono buone ragioni per affermare che a quel tempo la concentrazione 
%   plasmatica non sia uniforme: la sostanza non si el ancora distribuita 
%   nell intero compartimento plasmatico (PROBLEMI DI MESCOLAMENTO)

% IL MODELLO COMPARTIMENTALE HA COME IPOTESI QUELLA DI OMOGENEITA', QUINDI
% AL MINUTO 1 IL COMPARTIMENTO NON SI COMPORTA ANCORA IN MANIERA OMOGENEA E
% SI DEVE "ADATTARE"

c4 = CONCENTRAZIONE - BASALE;  % CONCENTRAZIONE - BASALE
C4 = c4(TEMPO>1);              % CONCENTRAZIONE - BASALE DAL MINUTO 2

% ESTRAGGO I DATI DAL FILE .dat

IMPORT = importdata('DatiCPsog4.dat');

dose = 49650;         % [pmol]    
CV = 0.04;            % CONSIDERO IL CASO A CV COSTANTE
T4 = TEMPO(TEMPO>1);  % TEMPO DAL MINUTO 2

figure(2)
subplot(2,2,1)
plot(TEMPO,CONCENTRAZIONE,'ko')
title('SCALA NORMALE')
xlabel('Tempo [minuti]')
ylabel('Concentrazione [pmol/ml]')

subplot(2,2,2)
semilogy(CONCENTRAZIONE)
title('originali scala semilog')

subplot(2,2,3)                  % ESCLUDO MINUTO 0 COL BASALE E MINUTO 1 
plot(T4,C4,'ko')
title('nuovi scala normale')
xlabel('Tempo [minuti]')
ylabel('Concentrazione [pmol/ml]')


subplot(2,2,4)
semilogy(C4)
title('nuovi scala semilog')


% MATRICE SIGMA (VARIANZA DELL'ERRORE DI MISURA)

% W PESO = ALL'INVERSO DELLA MATRICE SIGMA

%        DA EDB:
%               CV=(SD)/E[x]
%               var=SD^2

sigma = (CV*diag(CONCENTRAZIONE(TEMPO>1)))^2;   

% SIGMA V MATRICE DI COVARIANZA DELL'ERRORE DI MISURA 
% MANTENGO IL BASALE PERCHE' L'ERRORE E' PROPORZIONALE ALLA MISURA

% STANDARDIZZO MOLTIPLICANDO PER LA CONCENTRAZIONE IPOTIZZANDO VALORE MEDIO
% NELLE MISURE
% E NE CALCOLO LA VARIANZA

W = inv(sigma);  

% SIGMA SARA' DIAGONALE PERCHE' IPOTIZZO MISURE SCORRELATE QUINDI I TERMINI
% FUORI DALLA DIAGONALE SARANNO PARI A 0 

% CONSIDERO MISURE DI CONCENTRAZIONE COL BASALE ESSENDO LA VARIANZA
% OPERATORE CHE TRASCURA EFFETTI ADDITIVI
% VAR[Y+K] = VAR[Y]   K COSTANTE


%% %%%%%%%%----------%%%%%%%%%% STIMA LS %%%%%%%%%%%----------%%%%%%%%%%%%%

% 1 ESPONENZIALE

% func = @(x) ANDAMENTO DELLA VARIABILE

% p1 PARAMETRI STIMATI A MANO (INIZIALIZZO L'ALGORITMO)
% p_1 = PARAMETRI STIMATI DALLA lsqnonlin


%       A = intercetta su scala semilog (Y)

%      alpha = pendenza su scala semilog (TRA DUE PUNTI TROVO m COEFF
%              ANGOLARE)


% DAL GRAFICO LA CONCENTRAZIONE SI ESAURISCE CIRCA A 120
% TEMPO PARI A 5*TAU


tau1 = 120/5;   % FACCIO RETTA TRA INIZIO E FINE PERCHE' PRESUPPOSTO UN 
                % MONO- ESPONENZIALE MENTRE PER I 2 E 3 MR LO PORTO DIETRO
                % 120 E' TEMPO FINALE
alpha = 1/tau1;

p1 = [11,alpha]';    
exp1 = @(x)((x(1).*exp(-x(2).*T4))-C4);

% X(1) PARAMETRO A 
% X(2) PARAMETRO alpha  c4 = z

% FUNZIONE DEI MINIMI QUADRATI

[p_1,WRSS1,res1,exitflag1,output1,lambda1,jac1] = lsqnonlin(exp1,p1);

% FUNZIONE CON PARAMETRI STIMATI DALLA lsqnonlin

exp_1 = (p_1(1).*exp(-p_1(2).*T4));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% modello LS con 2 exp

%      A  alpha  B   beta  
p2 = [10, alpha, 1, 0.125]';   

exp2 = @(x)((x(1).*exp(-x(2).*T4))+(x(3).*exp(-x(4).*T4)))-C4;

% FUNZIONE DEI MINIMI QUADRATI

[p_2,WRSS2,res2,exitflag2,output2,lambda2,jac2] = lsqnonlin(exp2,p2);

% FUNZIONE CON PARAMETRI STIMATI DALLA lsqnonlin

exp_2 = ((p_2(1).*exp(-p_2(2).*T4))+(p_2(3).*exp(-p_2(4).*T4)));


% ---------------------------------------------------------------------------- %
% modello LS con 3 exp

%      A   alpha  B   beta    C   gamma
p3 = [9.5, alpha, 1, 0.0125, 0.5, 0.05]'; 

exp3 = @(x)((x(1).*exp(-x(2).*T4))+(x(3).*exp(-x(4).*T4))+(x(5).*exp(-x(6).*T4)))-C4;

% FUNZIONE DEI MINIMI QUADRATI

[p_3,WRSS3,res3,exitflag3,output3,lambda3,jac3] = lsqnonlin(exp3,p3);

% FUNZIONE CON PARAMETRI STIMATI DALLA lsqnonlin

exp_3 = (p_3(1).*exp(-p_3(2).*T4))+(p_3(3).*exp(-p_3(4).*T4))+(p_3(5).*exp(-p_3(6).*T4));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% GRAFICO I FITTING LS SIA SU SCALA NORMALE SIA SU SEMILOG

figure('Name','Modello LS')

subplot(3,2,1)                                          % 1 EXP
plot(T4,C4,'ko',T4,p_1(1).*exp(-p_1(2).*T4),'b-')
title('fit con 1 exp LS')
legend('Data','Fit')
xlabel('tempo [minuti]')
ylabel('concentrazioni [pmol/ml]')

subplot(3,2,2)                                          % 1 EXP SEMILOG
semilogy(T4,C4,'r',T4,p_1(1).*exp(-p_1(2).*T4),'b')
title('semilogaritmica con 1 exp LS')
legend('Data','Fit')
xlabel('tempo')
ylabel('concentrazioni scala log')



subplot(3,2,3)                                          % 2 EXP 
plot(T4,C4,'ko',T4,(p_2(1).*exp(-p_2(2).*T4))+(p_2(3).*exp(-p_2(4).*T4)),'b-')
title('fit con 2 exp LS')
legend('Data','Fit')
xlabel('tempo [minuti]')
ylabel('concentrazioni [pmol/ml]')

subplot(3,2,4)                                          % 2 EXP SEMILOG
semilogy(T4,C4,'r',T4,(p_2(1).*exp(-p_2(2).*T4))+(p_2(3).*exp(-p_2(4).*T4)),'b')
title('semilogaritmica con 2 exp LS')
legend('Data','Fit')
xlabel('tempo')
ylabel('concentrazioni scala log')



subplot(3,2,5)                                          % 3 EXP 
plot(T4,C4,'ko')
hold on
plot(T4,(p_3(1).*exp(-p_3(2).*T4))+(p_3(3).*exp(-p_3(4).*T4))+(p_3(5).*exp(-p_3(6).*T4)),'b-')
title('fit con 3 exp LS')
legend('Data','Fit')
xlabel('tempo [minuti]')
ylabel('concentrazioni [pmol/ml]')

subplot(3,2,6)                                          % 3 EXP SEMILOG
semilogy(T4,C4,'m')
hold on
semilogy(T4,(p_3(1).*exp(-p_3(2).*T4))+(p_3(3).*exp(-p_3(4).*T4))+(p_3(5).*exp(-p_3(6).*T4)),'b')
title('semilogaritmica con 3 exp LS')
legend('Data','Fit')
xlabel('tempo')
ylabel('concentrazioni scala log')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% ANALISI dei RESIDUI LS 

% VALUTO LA BIANCHEZZA CON 'ISPEZIONE GRAFICA'
% HO OVERFITTING SE I RESIDUI SONO TROPPO PICCOLI
% CONSIDERO BUONI I FITTING CON RESIDUI 
% COMPRESI AL 60% TRA -1 E 1
%          AL 90% TRA -2 E 2


% DAL GRAFICO DEI RESIDUI SU SCALA TEMPORALE VEDO CHE ESISTE UN TREND:
% PRIMA SONO TUTTI POSITIVI E POI TUTTI NEGATIVI

% E SOPRATTUTTO DA QUELLO CON X=CONCENTR Y=RESIDUI VEDO CHE ESISTE UNA
% CHIARA DIPENDENZA DALLA GRANDEZZA DELLA CONCENTRAZIONE: AL CRESCERE DELLA
% CONCENTRAZIONE CRESCONO ANCHE I RESIDUI

% IL TUTTO NON È COMPATIBILE CON LA HP CHE HO FATTO: 
% IL MODELLO CORRETTO NON È A SD COSTANTE

figure('Name','Residui caso LS')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

subplot(3,2,1)                                  % RESIDUI NEL TEMPO 1 EXP
plot(T4,res1)
title('residuo 1 exp LS')
axis([0 180 -2 2])
hold on
scatter(T4,res1,'ko','b')
hold on
yline(0)
yline(-1,'r--')
yline(1,'r--')
yline(-2,'r--')
yline(2,'r--')
legend('Interpolazione','Residui')
xlabel('tempo [minuti]')
ylabel('residui')

subplot(3,2,2)                                  % RES SU CONC 1 EXP
scatter(C4,res1,'ko','b')
title('RESIDUI IN FUNZIONE DI C 1 exp LS')
axis([0 15 -2 2])
hold on
yline(0)
xlabel('concentrazione [pmol/ml]')
ylabel('residui')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

subplot(3,2,3)                                  % RESIDUI NEL TEMPO 2 EXP
plot(T4,res2)
title('residuo 2 exp LS')
axis([0 180 -2 2])
hold on
scatter(T4,res2,'ko','b')
hold on
yline(0)
yline(-1,'r--')
yline(1,'r--')
yline(-2,'r--')
yline(2,'r--')
legend('Interpolazione','Residui')
xlabel('tempo [minuti]')
ylabel('residui')

subplot(3,2,4)                                  % RES SU CONC 2 EXP
scatter(C4,res2,'ko','b')
title('RESIDUI IN FUNZIONE DI C 2 exp LS')
axis([0 15 -2 2])
hold on
yline(0)
xlabel('concentrazione [pmol/ml]')
ylabel('residui')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

subplot(3,2,5)                                  % RESIDUI NEL TEMPO 3 EXP
plot(T4,res3)
title('residuo 3 exp LS')
axis([0 180 -2 2])
hold on
scatter(T4,res3,'ko','b')
hold on
yline(0)
yline(-1,'r--')
yline(1,'r--')
yline(-2,'r--')
yline(2,'r--')
legend('Interpolazione','Residui')
xlabel('tempo [minuti]')
ylabel('residui')

subplot(3,2,6)                                  % RES SU CONC 3 EXP
scatter(C4,res3,'ko','b')
title('RESIDUI IN FUNZIONE DI C 3 exp LS')
axis([0 15 -2 2])
hold on
yline(0)
xlabel('concentrazione [pmol/ml]')
ylabel('residui')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% INDICI DI MERITO: QUANTO STO FITTANDO BENE?

% N = numero dei dati (30)

% M = NUMERO DI PARAMETRI DA STIMARE

N = length(T4);

M1 = 2;     % p = [A alpha]
M2 = 4;     % p = [A alpha B beta]
M3 = 6;     % p = [A alpha B C gamma]

% CONSIDERO WRSS ESTRATTI DALLE lsqnolin

WRSS = [WRSS1 WRSS2 WRSS3]; 

M = [M1 M2 M3];

% INSERISCO LA MATRICE DI FISHER NEI 3 CASI PARI A Ht * sigmav^-1 *H
% sigma^-1 = W MA NEI MINIMI QUADRATI NON PESATI PONGO = 1

%---------------------------- 1 EXP ---------------------------------------
F1_LS = jac1'*jac1;

% MEDIA ARITMETICA DEGLI AUTOVALORI = TRACCIA FRATTO LA DIMENSIONE DELLA
%                                     MATRICE DI FISHER

la1LS = trace(F1_LS)/M1;

% MEDIA GEOMETRICA DEGLI AUTOVALORI = RADICE M-ESIMA DEL DETERMINANTE 

lg1LS = (det(F1_LS))^(1/M1);

% CALCOLO INDICE C1

c11LS = M1/(2*log(la1LS/lg1LS));

%---------------------------- 2 EXP ---------------------------------------


F2_LS = jac2'*jac2;
la2LS = trace(F2_LS)/M2;
lg2LS = (det(F2_LS))^(1/M2);
c12LS = M2/(2*log(la2LS/lg2LS));

%---------------------------- 3 EXP ---------------------------------------

F3_LS = jac3'*jac3;
la3LS = trace(F3_LS)/M3;
lg3LS = (det(F3_LS))^(1/M3);
c13LS = M3/(2*log(la3LS/lg3LS));


c1LS = [c11LS c12LS c13LS];


for i=1:3

    AIC(i) = WRSS(i)+2*M(i);  % CRITERIO AKAIKE
                              % MIGLIORE AL DIMINUIRE DI AIC 
                              % USATO CON N PICCOLO
 
    BIC(i) = WRSS(i)+log(N)*M(i);  % BAYESIAN INFORMTION CRITERION
                                   % DETTO ANCHE CRITERIO DI SCHWARTZ
                                   % USATO PER N GRANDI

   
   % MDL e FPE COMMENTATI(%) PERCHE' FUNZIONANTI CON TANTI DATI
   %    QUINDI NON NEI NOSTRI CASI DI FARMACOCINETICA
   
   % MDL(i) = WRSS(i)/N+log(N)*(M(i)/N); 
   % Minimum Description Length
   % FPE(i) = WRRS(i)(N+M(i))/(N-M(i))
   % Final Predictor Error

   ICOMP_LS(i) = WRSS(i)+(2*c1LS(i));  

end 


% NB. SE FACESSI GLI INDICI DI MERITO SU TUTTI I SOGGETTI
%     VEDO CHE ALCUNI RISULTANO ESSERE MODELLIZZATI MEGLIO
%     CON 2 EXP (ES. SOGGETTO 1,5,7) 
%     ALTRI CON 3 EXP (ES. SOGGETTO 2,3,4,6)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% AFFIDABILITA' DELLE STIME NEL CASO DI LS
% ANALIZZO LA MATRICE SIGMA P (MATRICE DI COVARIANZA DELLE STIME)

% IN QUESTO CASO NON LINEARE NEI PARAMETRI QUINDI LA lsqnolin RESTITUISCE
% LA MATRICE DI SENSITIVITA' S (gradiente)


% DA EDB  CV=SD/E[x]
%         var=SD^2
% VARIANZA SULLA DIAGONALE DELLA SIGMA P

% RICORDANDO CHE F^-1 = SIGMAP

% CV 1 exp LS
sigma_p1LS = (F1_LS)^-1;
diag_sigmap1LS = sigma_p1LS*eye(M1); % REALIZZA LA MATRICE DAGLI SPARSE DOUBLE
varp1LS = diag(diag_sigmap1LS);      % RITAGLIA GLI ELEMENTI SULLA DIAGONALE
CV_p1LS = (sqrt(varp1LS)'./abs(p_1))*100;

% CV 2 exp LS
sigma_p2LS = (F2_LS)^-1;
diag_sigmap2LS = sigma_p2LS*eye(M2);
varp2LS = diag(diag_sigmap2LS);
CV_p2LS = (sqrt(varp2LS)'./abs(p_2))*100;

% CV 3 exp LS
sigma_p3LS = (F3_LS)^-1;
diag_sigmap3LS = sigma_p3LS*eye(M3);
varp3LS = diag(diag_sigmap3LS);
CV_p3LS = (sqrt(varp3LS)'./abs(p_3))*100;

% MI ACCORGO CHE PER 1 E 2 ESPONENZIALI LA MATRICE E' PIU' O MENO BUONA
% TUTT'ALTRO PER 3 ESPONENZIALIM --> VALORI TROPPO GRANDI SU DIAG E FUORI


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%% %%%%%%%%----------%%%%%%%%%% STIMA WLS %%%%%%%%%%%----------%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% W PARI A sigma_v^-1
% VETTORE DI INIZIALIZZAZIONE = A LS

%---------------------------- WLS 1 EXP -----------------------------------

% UTILIZZIAMO LA RADICE DI W PERCHE' LA lsqnolin ELEVA L'ERRORE, QUINDI
% QUELLO CHE VOGLIAMO OTTENERE NOI E' W*(z-y)^2

exppes1 = @(x) (sqrt(W)*((x(1).*exp(-x(2).*T4))-C4));
[p_pes1,WRSSpes1,res_pes1,exitflag_pes1,output_pes1,lambda_pes1,jac_pes1]=lsqnonlin(exppes1,p1);
exppes_1 = (p_pes1(1).*exp(-p_pes1(2).*T4));


%---------------------------- WLS 2 EXP -----------------------------------


exppes2 = @(x) sqrt(W)*(((x(1).*exp(-x(2).*T4))+(x(3).*exp(-x(4).*T4))-C4));
[p_pes2,WRSSpes2,res_pes2,exitflag_pes2,output_pes2,lambda_pes2,jac_pes2]=lsqnonlin(exppes2,p2);
exppes_2 = ((p_pes2(1).*exp(-p_pes2(2).*T4))+(p_pes2(3).*exp(-p_pes2(4).*T4)));

%---------------------------- WLS 3 EXP -----------------------------------


exppes3 = @(x) sqrt(W)*((x(1).*exp(-x(2).*T4))+(x(3).*exp(-x(4).*T4))+(x(5).*exp(-x(6).*T4))-C4);
[p_pes3,WRSSpes3,res_pes3,exitflag_pes3,output_pes3,lambda_pes3,jac_pes3]=lsqnonlin(exppes3,p3);
exppes_3 = (p_pes3(1).*exp(-p_pes3(2).*T4))+(p_pes3(3).*exp(-p_pes3(4).*T4))+(p_pes3(5).*exp(-p_pes3(6).*T4));



figure('Name','Modello WLS')

subplot(3,2,1)                                          % 1 EXP WLS
plot(T4,C4,'ko',T4,p_pes1(1).*exp(-p_pes1(2).*T4),'b-')
title('fit con 1 exp WLS')
legend('Data','Fit')
xlabel('tempo [minuti]')
ylabel('concentrazioni [pmol/ml]')

subplot(3,2,2)                                          % 1 EXP WLS SEMILOG
semilogy(T4,C4,'r',T4,p_pes1(1).*exp(-p_pes1(2).*T4),'b')
title('semilogaritmica con 1 exp WLS')
legend('Data','Fit')
xlabel('tempo')
ylabel('concentrazioni scala log')



subplot(3,2,3)                                          % 2 EXP WLS
plot(T4,C4,'ko',T4,(p_pes2(1).*exp(-p_pes2(2).*T4))+(p_pes2(3).*exp(-p_pes2(4).*T4)),'b-')
title('fit con 2 exp WLS')
legend('Data','Fit')
xlabel('tempo [minuti]')
ylabel('concentrazioni [pmol/ml]')

subplot(3,2,4)                                          % 2 EXP WLS SEMILOG
semilogy(T4,C4,'r',T4,(p_pes2(1).*exp(-p_pes2(2).*T4))+(p_pes2(3).*exp(-p_pes2(4).*T4)),'b')
title('semilogaritmica con 2 exp WLS')
legend('Data','Fit')
xlabel('tempo')
ylabel('concentrazioni scala log')



subplot(3,2,5)                                          % 3 EXP WLS
plot(T4,C4,'ko')
hold on
plot(T4,(p_pes3(1).*exp(-p_pes3(2).*T4))+(p_pes3(3).*exp(-p_pes3(4).*T4))+(p_pes3(5).*exp(-p_pes3(6).*T4)),'b-')
title('fit con 3 exp WLS')
legend('Data','Fit')
xlabel('tempo [minuti]')
ylabel('concentrazioni [pmol/ml]')

subplot(3,2,6)                                          % 3 EXP WLS SEMILOG
semilogy(T4,C4,'m')
hold on
semilogy(T4,(p_pes3(1).*exp(-p_pes3(2).*T4))+(p_pes3(3).*exp(-p_pes3(4).*T4))+(p_pes3(5).*exp(-p_pes3(6).*T4)),'b')
title('semilogaritmica con 3 exp WLS')
legend('Data','Fit')
xlabel('tempo')
ylabel('concentrazioni scala log')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% ANALISI dei RESIDUI WLS

% VALUTO LA BIANCHEZZA CON 'ISPEZIONE GRAFICA'
% HO OVERFITTING SE I RESIDUI SONO TROPPO PICCOLI
% CONSIDERO BUONI I FITTING CON RESIDUI 
% COMPRESI AL 60% TRA -1 E 1
%          AL 90% TRA -2 E 2

% I RESIDUI VANNO 'UN PO' SU UN PO' GIÙ'
% HANNO CIRCA LA STESSA AMPIEZZA E NON VEDO SISTEMATICITÀ

% IL MODELLO DELL'ERRORE CON CV COSTANTE SEMBRA RAGIONEVOLE



figure('Name','Residui caso WLS')

subplot(3,2,1)                           % RESIDUI WLS NEL TEMPO 1 EXP
plot(T4,res_pes1)
title('residuo norm 1 exp WLS')
hold on
scatter(T4,res_pes1,'ko','b')
hold on
yline(0)
yline(-1,'r--')
yline(1,'r--')
yline(-2,'r--')
yline(2,'r--')
legend('Interpolazione','Residui normalizzati')
xlabel('tempo [minuti]')
ylabel('residui norm')


subplot(3,2,2)                           % RES SU CONC WLS 1 EXP
scatter(C4,res_pes1,'ko','b')
title('RESIDUI IN FUNZIONE DI C norm 1 exp WLS')
axis([0 15 -2 2])
hold on
yline(0)
xlabel('concentrazione [pmol/ml]')
ylabel('residui norm')



subplot(3,2,3)                           % RESIDUI WLS NEL TEMPO 2 EXP
plot(T4,res_pes2)
title('residuo norm 2 exp WLS')
hold on
scatter(T4,res_pes2,'ko','b')
hold on
yline(0)
yline(-1,'r--')
yline(1,'r--')
yline(-2,'r--')
yline(2,'r--')
legend('Interpolazione','Residui normalizzati')
xlabel('tempo [minuti]')
ylabel('residui norm')


subplot(3,2,4)                           % RES SU CONC WLS 2 EXP
scatter(C4,res_pes2,'ko','b')
title('RESIDUI IN FUNZIONE DI C norm 2 exp WLS')
hold on
yline(0)
xlabel('concentrazione [pmol/ml]')
ylabel('residui norm')



subplot(3,2,5)                           % RESIDUI WLS NEL TEMPO 3 EXP
plot(T4,res_pes3)
title('residuo norm 3 exp WLS')
hold on
scatter(T4,res_pes3,'ko','b')
hold on
yline(0)
yline(-1,'r--')
yline(1,'r--')
yline(-2,'r--')
yline(2,'r--')
legend('Interpolazione','Residui normalizzati')
xlabel('tempo [minuti]')
ylabel('residui')


subplot(3,2,6)                           % RES SU CONC WLS 3 EXP
scatter(C4,res_pes3,'ko','b')
title('RESIDUI IN FUNZIONE DI C norm 3 exp WLS')
hold on
yline(0)
xlabel('concentrazione [pmol/ml]')
ylabel('residui norm')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% INDICI DI MERITO: QUANTO STO FITTANDO BENE?

% N = numero dei dati (30)

% M = NUMERO DI PARAMETRI DA STIMARE


WRSSpes = [WRSSpes1 WRSSpes2 WRSSpes3]'; % IN OUTPUT DALLA lsqnolin


F1 = jac_pes1'*W*jac_pes1;
la1 = trace(F1)/M1;
lg1 = (det(F1))^(1/M1);
c11 = M1/(2*log(la1/lg1));

F2 = jac_pes2'*W*jac_pes2;
la2 = trace(F2)/M2;
lg2 = (det(F2))^(1/M2);
c12 = M2/(2*log(la2/lg2));

F3 = jac_pes3'*W*jac_pes3;
la3 = trace(F3)/M3;
lg3 = (det(F3))^(1/M3);
c13 = M3/(2*log(la3/lg3));

c1 = [c11 c12 c13];

for i=1:3

    AICpes(i) = WRSSpes(i)+2*M(i);  
 
    BICpes(i) = WRSSpes(i)+log(N)*M(i);  

    ICOMP(i) = WRSSpes(i)+(2*c1(i));  
   
  
   
   % MDLpes(i) = WRSSpes(i)/N+log(N)*(M(i)/N); 
   % Minimum Description Length
   % FPEpes(i) = WRRSpes(i)(N+M(i))/(N-M(i))
   % Final Predictor Error

end 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% VALUTAZIONI GRAFICHE confronto STIME LS e WLS

% COME SI DISTRIBUISCONO I DATI E LE PREDIZIONI = COVARIANZA
% ASSE X: DATI
%      Y: PREDIZIONI


figure(7)
subplot(3,2,1)                                          % 1 EXP LS COV
scatter(C4,p_1(1).*exp(-p_1(2).*T4),'ko','b')
hold on
plot(TEMPO,TEMPO)
axis([0 12 0 12])
title('1 exp LS')
xlabel('dati')
ylabel('predizioni')


subplot(3,2,3)                                         % 2 EXP LS COV

scatter(C4,(p_2(1).*exp(-p_2(2).*T4))+(p_2(3).*exp(-p_2(4).*T4)),'ko','b')
hold on
plot(TEMPO,TEMPO)
axis([0 12 0 15])
title('2 exp LS')
xlabel('dati')
ylabel('predizioni')


subplot(3,2,5)                                         % 3 EXP LS COV

scatter(C4,(p_3(1).*exp(-p_3(2).*T4))+(p_3(3).*exp(-p_3(4).*T4))+(p_3(5).*exp(-p_3(6).*T4)),'ko','b')
hold on
plot(TEMPO,TEMPO)
axis([0 12 0 12])
title('3 exp LS')
xlabel('dati')
ylabel('predizioni')



subplot(3,2,2)                                         % 1 EXP WLS COV
scatter(C4,p_pes1(1).*exp(-p_pes1(2).*T4),'ko','b')
hold on
plot(TEMPO,TEMPO)
axis([0 12 0 7])
title('1 exp WLS')
xlabel('dati')
ylabel('predizioni')


subplot(3,2,4)                                         % 2 EXP WLS COV

scatter(C4,(p_pes2(1).*exp(-p_pes2(2).*T4))+(p_pes2(3).*exp(-p_pes2(4).*T4)),'ko','b')
hold on
plot(TEMPO,TEMPO)
axis([0 12 0 12])
title('2 exp WLS')
xlabel('dati')
ylabel('predizioni')


subplot(3,2,6)                                         % 3 EXP WLS COV

scatter(C4,(p_pes3(1).*exp(-p_pes3(2).*T4))+(p_pes3(3).*exp(-p_pes3(4).*T4))+(p_pes3(5).*exp(-p_pes3(6).*T4)),'ko','b')
hold on
plot(TEMPO,TEMPO)
axis([0 12 0 15])
title('3 exp WLS')
xlabel('dati')
ylabel('predizioni')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% HISTOGRAMMI e GAUSSIANA

%PLOTTO LA GAUSSIANA PER CONFRONTARE I RESIDUI

mu = 0;
sig = 1;
x_g = -15:0.1:15;
Gauss = exp(-(x_g-mu).^2./(2*sig^2))./(sig*sqrt(2*pi));



figure('Name','hist/gauss RESIDUI')

subplot(3,2,1)
histogram(res1,10,'Normalization','pdf')
axis([-3 3 0 0.8])
hold on
plot(x_g,Gauss,'r')
title('hist/gauss 1 exp LS')
legend('hist','gauss')


subplot(3,2,3)
histogram(res2,10,'Normalization','pdf')
axis([-4 2.5 0 1.5])
hold on
plot(x_g,Gauss,'r')
title('hist/gauss 2 exp LS')
legend('hist','gauss')


subplot(3,2,5)
histogram(res3,10,'Normalization','pdf')
axis([-4 2 0 1.5])
hold on
plot(x_g,Gauss,'r')
title('hist/gauss 3 exp LS')
legend('hist','gauss')


subplot(3,2,2)
histogram(res_pes1,10,'Normalization','pdf')
hold on
plot(x_g,Gauss,'r')
title('hist/gauss 1 exp WLS norm')
legend('hist','gauss')


subplot(3,2,4)
histogram(res_pes2,10,'Normalization','pdf')
hold on
plot(x_g,Gauss,'r')
title('hist/gauss 2 exp WLS norm')
legend('hist','gauss')


subplot(3,2,6)
histogram(res_pes3,10,'Normalization','pdf')
hold on
plot(x_g,Gauss,'r')
title('hist/gauss 3 exp WLS norm')
legend('hist','gauss')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% AFFIDABILITA' delle stime caso WLS

% CALCOLO IL NUOVO COEFFICIENTE DI VARIAZIONE CHE MI FORNISCE UNA MISURA
% DELL'AFFIDABILITÀ DELLE STIME:
% LA VARIANZA È DATA DAGLI ELEMENTI SULLA DIAGONALE DI SIGMA P

% CV 1 exp WLS
sigma_p1 = F1^(-1);
diag_sigmap1 = sigma_p1*eye(M1);
varp1 = diag(diag_sigmap1);
CV_p1 = (sqrt(varp1)'./abs(p_pes1))*100;

% CV 2 exp WLS
sigma_p2 = [F2]^(-1);
diag_sigmap2 = sigma_p2*eye(M2);
varp2 = diag(diag_sigmap2);
CV_p2 = (sqrt(varp2)'./abs(p_pes2))*100;

% CV 3 exp WLS
sigma_p3 = [F3]^(-1);
diag_sigmap3 = sigma_p3*eye(M3);
varp3 = diag(diag_sigmap3);
CV_p3 = (sqrt(varp3)'./abs(p_pes3))*100;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% vedo che la deviazione standard dei residui non è uguale a 1
% quindi ipotizzo che il CV non sia 0,04 ma magari 0,06

std1 = std(res_pes1);
std2 = std(res_pes2);
std3 = std(res_pes3);

% STIMO NUOVO CV TRAMITE ALPHA
% RMB. alpha = WRSS/(N-M)    N-M GRADI DI LIBERTA'
%      alpha = (CV)^2       

for i=1:3
    xci(i) = WRSSpes(i)/(N-M(i));
end 

for i=1:3
    CV_new(i) = sqrt(xci(i));
end


%% %%%%%%%%%%%%%%%%%%%%%% MODELLI DI POPOLAZIONE %%%%%%%%%%%%%%%%%%%%%%%% %

% CARICO TUTTI I DATI DI TUTTI I SOGGETTI RIMUOVENDO IL BASALE E RIMUOVENDO
% MISURE AL MINUTO 0 E 1

TEMPO = S1(:,1);         
conc1 = S1(:,2);         
basale1 = conc1(TEMPO==0);
cc1 = conc1(TEMPO>1);
c1 = conc1 - basale1;
C1 = c1(TEMPO>1);

conc2 = S2(:,2);         
basale2 = conc2(TEMPO==0);
c2 = conc2 - basale2;
cc2 = conc2(TEMPO>1);
C2 = c2(TEMPO>1);

conc3 = S3(:,2);         
basale3 = conc3(TEMPO==0);
c3 = conc3 - basale3;
cc3 = conc3(TEMPO>1);
C3 = c3(TEMPO>1);

conc4 = S4(:,2);         
basale4 = conc4(TEMPO==0);
c4 = conc4 - basale4;
cc4 = conc4(TEMPO>1);
C4 = c4(TEMPO>1);

conc5 = S5(:,2);         
basale5 = conc5(TEMPO==0);
c5 = conc5 - basale5;
cc5 = conc5(TEMPO>1);
C5 = c5(TEMPO>1);
     
conc6 = S6(:,2);         
basale6 = conc6(TEMPO==0);
c6 = conc6 - basale6;
cc6 = conc6(TEMPO>1);
C6 = c6(TEMPO>1);
  
conc7 = S7(:,2);         
basale7 = conc7(TEMPO==0);
c7 = conc7 - basale7;
cc7 = conc7(TEMPO>1);
C7 = c7(TEMPO>1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%% METODO NAIVE AVERAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%% %
       %     CONSIDERO SOLO IL MODELLO WLS CON 2 EXP


concAV = [C1 C2 C3 C4 C5 C6 C7]';  %% SENZA BASALE PER CALCOLARMI LA MEDIA
av = mean(concAV)';

wAV = [cc1 cc2 cc3 cc4 cc5 cc6 cc7]'; %% CON BASALE PER CALCOLARMI LA W
wav = mean(wAV)';
sigmaAV = (CV*diag(wav))^2;
WAV = inv(sigmaAV);         

NaiveAV = @(x) sqrt(WAV)*(((x(1).*exp(-x(2).*T4))+(x(3).*exp(-x(4).*T4))-av));
[p_pesAV,WRSSpesAV,res_pesAV,exitflag_pesAV,output_pesAV,lambda_pesAV,jac_pesAV]=lsqnonlin(NaiveAV,p2);
Naive_AV = ((p_pesAV(1).*exp(-p_pesAV(2).*T4))+(p_pesAV(3).*exp(-p_pesAV(4).*T4)));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%% METODO NAIVE POOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
       %     CONSIDERO SOLO IL MODELLO WLS CON 2 EXP

POOL = zeros(210,2);

for i=1:length(T4)

    if i==1

        pool = [T4(i) C1(i);
                T4(i) C2(i);
                T4(i) C3(i);
                T4(i) C4(i);
                T4(i) C5(i);
                T4(i) C6(i);
                T4(i) C7(i)];
   
        POOL = [pool];

    else

        pool = [T4(i) C1(i);
            T4(i) C2(i);
            T4(i) C3(i);
            T4(i) C4(i);
            T4(i) C5(i);
            T4(i) C6(i);
            T4(i) C7(i)];
   
        POOL = [POOL; pool];
    end
    

end 

TP = POOL(:,1);
CP = POOL(:,2);
sigmaPOOL = (CV*diag(CP))^2;    
WPOOL = inv(sigmaPOOL);     

NaivePOOL = @(x) sqrt(WPOOL)*(((x(1).*exp(-x(2).*TP))+(x(3).*exp(-x(4).*TP))-CP));
[p_pesPOOL,WRSSpesPOOL,res_pesPOOL,exitflag_pesPOOL,output_pesPOOL,lambda_pesPOOL,jac_pesPOOL]=lsqnonlin(NaivePOOL,p2);
Naive_POOL = ((p_pesPOOL(1).*exp(-p_pesPOOL(2).*TP))+(p_pesPOOL(3).*exp(-p_pesPOOL(4).*TP)));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%% METODO 2 STAGES STANDARD %%%%%%%%%%%%%%%%%%%%%%%% %

% mean(parametri) = i parametri di popolazione sono calcolati come media dei parametri
%                   stimati sull'intera popolazione
% cov(parametri) = considero qua anche la variabilità interindividuale
clear TwoSS

for i=1:7

    sigma2SS = (CV*diag(wAV(i,:)))^2;    
    W2SS = inv(sigma2SS);    

    % come partenza dell'algoritmo gli do Naive Average
    TwoStageStand = @(x) sqrt(W2SS)*(((x(1).*exp(-x(2).*T4))+(x(3).*exp(-x(4).*T4))-concAV(i,:)));
    [p_pesSS,WRSSpesSS,res_pesSS,exitflag_pesSS,output_pesSS,lambda_pesSS,jac_pesSS]=lsqnonlin(TwoStageStand,p_pesAV,[0 0 0 0]);
    
    TwoSS(:,i) = [p_pesSS(1) p_pesSS(2) p_pesSS(3) p_pesSS(4)]; 
end

TwoSS = abs(TwoSS)';
parTwoSS = mean(TwoSS);
parTwoSS = parTwoSS';
Two_SS = ((parTwoSS(1).*exp(-parTwoSS(2).*T4))+(parTwoSS(3).*exp(-parTwoSS(4).*T4)));


var_interind = cov(TwoSS);

%% PLOT DEI GRAFICI di POPOLAZIONE

figure('Name','Modelli di Popolazione')

%%%%%%%%%%%% Naive Average
subplot(3,2,1)
plot(T4,av,'*','Color',[0.9290 0.6940 0.1250])
hold on
plot(T4,Naive_AV,'b')
hold on
plot(T4,C4,'o')
title('Naive Average')
legend('dati medi','fit','pz.4')

subplot(3,2,2)
semilogy(av,'*','Color',[0.9290 0.6940 0.1250])
hold on
semilogy(Naive_AV,'b')
hold on
semilogy(C4,'o')
title('Naive Average Semilog')
legend('dati medi','fit','pz.4')

%%%%%%%%%%%% POOL
subplot(3,2,3)
plot(TP,CP,'*','Color',[0.9290 0.6940 0.1250])
hold on
plot(TP,Naive_POOL,'b')
title('Naive Poolled')
legend('dati','fit')

subplot(3,2,4)
semilogy(CP,'*','Color',[0.9290 0.6940 0.1250])
hold on
semilogy(Naive_POOL,'b')
title('Naive Poolled')
legend('dati','fit')

%%%%%%%%%%%% 2-stage standard
subplot(3,2,5)
plot(TP,CP,'*','Color',[0.9290 0.6940 0.1250])
hold on
plot(T4,Two_SS,'b')
title('2 Stage Standard')
legend('dati','fit')

subplot(3,2,6)
semilogy(CP,'*','Color',[0.9290 0.6940 0.1250])
hold on
semilogy(Two_SS,'b')
title('2 Stage Standard Semilog')
legend('dati','fit','pz.4')



%% REPORT CON DATI SIGNIFICATIVI
EXP = [1 2 3];

fid = fopen('ES3.txt','w');

fprintf(fid,'\n DATI RILEVANTI ESERCITAZIONE 3 \n');

fprintf(fid,'\n modello LS \n');
fprintf(fid,'\nexp    A        alpha          B             beta             C             gamma');
fprintf(fid,'\n%d\t%f\t%f\t%d\t%d\t%d\t%d',1,p_1(1),p_1(2));
fprintf(fid,'\n%d\t%f\t%f\t%d\t%d\t%d\t%d',2,p_2(3),p_2(4),p_2(1),p_2(2));
fprintf(fid,'\n%d\t%f\t%f\t%d\t%d\t%d\t%d',3,p_3(1),p_3(2),p_3(3),p_3(4),p_3(5),p_3(6));

fprintf(fid,'\n\n modello WLS \n');
fprintf(fid,'\nexp    A        alpha          B             beta             C             gamma');
fprintf(fid,'\n%d\t%f\t%f\t%d\t%d\t%d\t%d',1,p_pes1(1),p_pes1(2));
fprintf(fid,'\n%d\t%f\t%f\t%d\t%d\t%d\t%d',2,p_pes2(3),p_pes2(4),p_pes2(1),p_pes2(2));
fprintf(fid,'\n%d\t%f\t%f\t%d\t%d\t%d\t%d',3,p_pes3(5),p_pes3(6),p_pes3(1),p_pes3(2),p_pes3(3),p_pes3(4));


fprintf(fid,'\n\n\n\n analisi del WRSS \n');
fprintf(fid,'\nexp      WRSS           WRSS pesato');
for i=1:3
    fprintf(fid, '\n%d\t\t%f\t\t%f',EXP(i),WRSS(i),WRSSpes(i));
     
end 


fprintf(fid,'\n\n\n\n analisi delle cifre di merito \n');
fprintf(fid,'\nexp      AIC             AIC pesato         BIC         BIC pesato          ICOMP         ICOMP pesato');
for i=1:3
    fprintf(fid, '\n%d\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f',EXP(i),AIC(i),AICpes(i),BIC(i),BICpes(i),ICOMP_LS(i),ICOMP(i));
     
end 


fprintf(fid,'\n\n\n\n nuovi CV delle STIME dei PARAMETRI \n');
fprintf(fid,'\nexp    modello        A          alpha            B             beta             C             gamma');
fprintf(fid, '\n%d\t\t%s\t\t%f\t\t%f\t',1,'LS',CV_p1LS(1,1),CV_p1LS(2,2));
fprintf(fid, '\n%d\t\t%s\t\t%f\t\t%f\t',1,'WLS',CV_p1(1,1),CV_p1(2,2));
fprintf(fid, '\n%d\t\t%s\t\t%f\t\t%f\t\t%f\t\t%f\t',2,'LS',CV_p2LS(1,1),CV_p2LS(2,2),CV_p2LS(3,3),CV_p2LS(4,4));
fprintf(fid, '\n%d\t\t%s\t\t%f\t\t%f\t\t%f\t\t%f\t',2,'WLS',CV_p2(1,1),CV_p2(2,2),CV_p2(3,3),CV_p2(4,4));
fprintf(fid, '\n%d\t\t%s\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t',3,'LS',CV_p3LS(1,1),CV_p3LS(2,2),CV_p3LS(3,3),CV_p3LS(4,4),CV_p3LS(5,5),CV_p3LS(6,6));
fprintf(fid, '\n%d\t\t%s\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t',3,'WLS',CV_p3(1,1),CV_p3(2,2),CV_p3(3,3),CV_p3(4,4),CV_p3(5,5),CV_p3(6,6));

fprintf(fid,'\n\n\n\n nuovi CV delle STIME caso WLS \n');
fprintf(fid,'\nexp      CV');
for i=1:3
    fprintf(fid, '\n%d\t\t%f\t',EXP(i),CV_new(i));
     
end 

fprintf(fid,'\n\n\n\n Modelli di Popolazione \n');
fprintf(fid,'\nmodello        A          alpha          B         beta');
fprintf(fid,'\n%s\t\t\t%f\t%f\t%f\t%f\t','AV',p_pesAV(3),p_pesAV(4),p_pesAV(1),p_pesAV(2));
fprintf(fid,'\n%s\t\t%f\t%f\t%f\t%f\t','Pooled',p_pesPOOL(3),p_pesPOOL(4),p_pesPOOL(1),p_pesPOOL(2));
fprintf(fid,'\n%s\t\t%f\t%f\t%f\t%f\t','WLS-2',parTwoSS(3),parTwoSS(4),parTwoSS(1),parTwoSS(2));


fprintf(fid,'\n\n\n\n\n');
fclose(fid);

% CV OTTIMI SE MINORI DI 10


