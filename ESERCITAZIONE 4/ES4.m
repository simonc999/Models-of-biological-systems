%% ESERCITAZIONE 4 MODELLI

% ATTENZIONE a fine script viene stmpato un REPORT riassuntivo
%            titolo del report = 'eserc4modelli.txt'
clc
clear 
close all
%% IDENTIFICAZIONE DI MODELLI COMPARTIMENTALI LINEARI


S4 = table2array(readtable('DatiCPsog4.dat'));

tempo = S4(:,1);         % [minuti]
conc = S4(:,2);          % [pmol/ml]
basale = conc(tempo==0);
c4 = conc - basale;
C4 = c4(tempo>1);
T4 = tempo(tempo>1);

CV = 0.04;               % caso a CV costante
sigma = (CV*diag(C4))^2;
W = inv(sigma);  
W_lsqnonlin = sqrt(inv(sigma));   

% NEL CASO NON PESATO INFATTI W È LA MATRICE IDENTITÀ
% AL POSTO DI RISCRIVERE LA FUNZIONE LR SENZA PESO, MI CONVIENE DIRLE CHE IL
% PESO È PARI A 1

W1 = eye(length(C4));

dose = 49650;     %[pmol]



%% da ESERCITAZIONE 1
% USO IL MODELLO BICOMPARTIMENTALE IMPLEMENTATO CON LA FUNZIONE LR
% (CON LR = LINEAR RESIDUAL) E LA LSQNONLIN PER RICAVARE I PARAMETRI
% INCOGNITI A PARTIRE DAI DATI 

% COME PARAMETRI INIZIALI DELL'ALGORITMO METTO QUELLI DELLA ESERCITAZIONE 1
V12 = 3.290;           %[litri]
k012 = 6.54e-2;        %[1/minuti]
k122 = 5.68e-2;        %[1/minuti]
k212 = 7.16e-2;        %[1/minuti]

p2comp = [k012,k122,k212,V12];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% modello LS

% LA FUNZIONE RESTITUISCE I RESIDUI IN FUNZIONE DEI PUNTATORI CHE SONO I
% PARAMETRI

fun2compLS = @(p2_compLS) LR(p2_compLS,T4,C4,W1,dose);
[p2_compLS,WRSS_2compLS,res_2compLS,exitflag_pes2cLS,output_pes2cLS,lambda_pes2cLS,jac_pes2cLS]=lsqnonlin(fun2compLS,p2comp);
[res2compLS,c2compLS,sys2compLS] = LR(p2_compLS,T4,C4,W1,dose);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% modello WLS

fun2comp = @(p2_comp) LR(p2_comp,T4,C4,W_lsqnonlin,dose);
[p2_comp,WRSS_2comp,res_2comp,exitflag_pes2c,output_pes2c,lambda_pes2c,jac_pes2c]=lsqnonlin(fun2comp,p2comp);
[res2comp,c2comp,sys2comp] = LR(p2_comp,T4,C4,W_lsqnonlin,dose);


% CV delle stime di K01 K12 K21 V 
M2 = 4;
F_2comp = jac_pes2c'*jac_pes2c;
sigma_p2comp = (F_2comp)^-1;
diag_sigmap2comp = sigma_p2comp*eye(M2);
varp2comp = diag(diag_sigmap2comp);
CV_2comp = (sqrt(varp2comp)'./abs(p2_comp))*100;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% GRAFICI MODELLO LS e WLS
% DEVO SCEGLIERE IL CANDIDATO CHE VADO POI A COMPARARE 
% CON IL MODELLO A 2 EXP DELLA ESERCITAZIONE 3


figure('Name','Scelta del modello 2-compartimenti')

subplot(2,3,1)                                          % 2 COMP LS
plot(T4,C4,'ko',T4,c2compLS,'b-')
title('fit LS')
legend('Data','Fit')
xlabel('tempo [min]')
ylabel('concentrazione [mg/l]')


subplot(2,3,2)                                          % 2 COMP LS SEMILOG
semilogy(T4,C4,'r',T4,c2compLS,'b')
title('semilog LS')
legend('Data','Fit')
xlabel('tempo [min]')
ylabel('concentrazione scala log [mg/l]')


subplot(2,3,3)                                          % 2 COMP LS COV
scatter(C4,c2compLS,'ko','b')
hold on
plot(tempo,tempo)
axis([0 12 0 12])
title('dati/predizioni LS')
xlabel('dati')
ylabel('predizioni')


subplot(2,3,4)                                          % 2 COMP WLS
plot(T4,C4,'ko',T4,c2comp,'b-')
title('fit WLS')
legend('Data','Fit')
xlabel('tempo [min]')
ylabel('concentrazione [mg/l]')


subplot(2,3,5)                                         % 2 COMP WLS SEMILOG
semilogy(T4,C4,'r',T4,c2comp,'b')
title('semilog WLS')
legend('Data','Fit')
xlabel('tempo [min]')
ylabel('concentrazione scala log [mg/l]')


subplot(2,3,6)                                          % 2 COMP WLS COV
scatter(C4,c2comp,'ko','b')
hold on
plot(tempo,tempo)
axis([0 12 0 12])
title('dati/predizioni WLS')
xlabel('dati')
ylabel('predizioni')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% ANALISI DEI RESIDUI LS E WLS
% EMERGE CHIARAMENTE CHE IL MODELLO A CV COSTANTE SIA IL MIGLIORE. 
% I RESIDUI NON PESATI INFATTI HANNO UN TREND DEFINITO: LA LORO AMPIEZZA
% CRESCE CON L'AUMENTARE DELLA CONCENTRAZIONE E NON HANNO DUNQUE
% DISTRIBUZIONE NORMALE (GAUSSIANA)

% UTILIZZO QUINDI IL METODO A CV COSTANTE

figure('Name','Analisi dei residui')



subplot(2,3,1)                                 % RESIDUI NEL TEMPO 2C LS
plot(T4,res2compLS)
title('residuo LS')
axis([0 180 -3 3])
hold on
scatter(T4,res2compLS,'ko','b')
hold on
yline(0)
yline(-1,'r--')
yline(1,'r--')
yline(-2,'r--')
yline(2,'r--')
legend('Interpolazione','Residui')
xlabel('tempo [minuti]')
ylabel('residui norm')


subplot(2,3,4)                                 % RESIDUI NEL TEMPO 2C WLS
plot(T4,res2comp)
title('residuo norm WLS')
hold on
scatter(T4,res2comp,'ko','b')
hold on
yline(0)
yline(-1,'r--')
yline(1,'r--')
yline(-2,'r--')
yline(2,'r--')
legend('Interpolazione','Residui normalizzati')
xlabel('tempo [minuti]')
ylabel('residui norm')




subplot(2,3,2)                                  % RES SU CONC 2C LS
scatter(C4,res2compLS,'ko','b')
title('RESIDUI IN FUNZIONE DI C LS')
axis([0 15 -2 2])
hold on
yline(0)
xlabel('concentrazione [pmol/ml]')
ylabel('residui')


subplot(2,3,5)                                  % RES SU CONC 2C WLS
scatter(C4,res2comp,'ko','b')
title('RESIDUI IN FUNZIONE DI C norm WLS')
hold on
yline(0)
xlabel('concentrazione [pmol/ml]')
ylabel('residui norm')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% HISTOGRAMMI e GAUSSIANA 


mu = 0;
sig = 1;
x_g = -15:0.1:15;
Gauss = exp(-(x_g-mu).^2./(2*sig^2))./(sig*sqrt(2*pi));

% LS
subplot(2,3,3)
histogram(res2compLS,10,'Normalization','pdf')
axis([-4 3 0 0.6])
hold on
plot(x_g,Gauss,'r')
title('hist/gauss LS')
legend('hist','gauss')

% WLS
subplot(2,3,6)
histogram(res2comp,10,'Normalization','pdf')
hold on
plot(x_g,Gauss,'r')
title('hist/gauss WLS')
legend('hist','gauss')
sgtitle('ANALISI RESIDUI CON PARAMETRI STIMATI DAL MODELLO COMPARTIMENTALE');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% RICAVO I PARAMETRI dal modello 2-comp WLS E CONVERTO ALLA ESPONENZIALE

% ESTRAGGO I PARAMETRI FORNITI IN OUTPUT DALLA lsqnolin CON IL METODO WLS

V_stim = p2_comp(4);
K01_stim = abs(p2_comp(1));
K12_stim = p2_comp(2);
K21_stim = p2_comp(3);

% RICAVO [A alpha B beta] 

transfunc = tf(sys2comp);
num = [14.03 0.7793];
den = [1 0.1833 0.003808];

% con 'residue' che converts the partial fraction expansion back to 
% the ratio of two polynomials and returns the coefficients in b and a.

[r,p] = residue(num,den);
Alpha_2comp = abs(p(1));
Beta_2comp = abs(p(2));
A_2comp = r(1);
B_2comp = r(2);



%% da ESERCITAZIONE 3
% RICAVO NUOVAMENTE I PARAMETRI DELLA ESERCITAZIONE 3 COSÌ DA CONFRONTARLI
% CON QUELLI OTTENUTI DAL MODELLO COMPARTIMENTALE

% MODELLO WLS CON 2 EXP
% CONSIDERO IL MODELLO A DUE ESPONENZIALI PERCHÈ SUGGERITO NEL TESTO
% CONSIDERO IL MODELLO PESATO PERCHÈ ABBIAMO VISTO NELLA ESERC 3 CHE IL
%           MODELLO DELL'ERRORE A CV COSTANTE È QUELLO PIÙ CORRETTO

tau1 = 120/5;
alpha = 1/tau1;

p2 = [10,alpha,1,0.125]';  
exppes2 = @(x) (W_lsqnonlin)*(((x(1).*exp(-x(2).*T4))+(x(3).*exp(-x(4).*T4))-C4));
[p_pes2,WRSSpes2,res_pes2,exitflag_pes2,output_pes2,lambda_pes2,jac_pes2]=lsqnonlin(exppes2,p2);
exppes_2 = ((p_pes2(1).*exp(-p_pes2(2).*T4))+(p_pes2(3).*exp(-p_pes2(4).*T4)));

% LA lsqnolin RESTITUIRA' I PARAMETRI DELL'ESPONENZIALE DOPPIO

A_stim = p_pes2(3);
Alpha_stim = p_pes2(4);
B_stim = p_pes2(1);
Beta_stim = p_pes2(2);


%% GRAFICI COMPARATIVI ESERC 1 e ESERC 3

%        2-COMP_WLS                 VS                 2 exp WLS

figure('Name','2-COMP e 2 exp WLS')

subplot(2,3,1)                                      % 2 COMP
plot(T4,C4,'ko',T4,c2comp,'b-')
title('fit con 2-compartimenti')
legend('Data','Fit')
xlabel('tempo [min]')
ylabel('concentrazione [mg/l]')


subplot(2,3,2)                                      % 2 COMP SEMILOG
semilogy(T4,C4,'r',T4,c2comp,'b')
title('semilogaritmica con 2-compartimenti')
legend('Data','Fit')
xlabel('tempo [min]')
ylabel('concentrazione scala log [mg/l]')


subplot(2,3,3)                                      % 2 COMP COV
scatter(C4,c2comp,'ko','b')
hold on
plot(tempo,tempo)
axis([0 12 0 12])
title('2-compartimenti')
xlabel('dati')
ylabel('predizioni')


subplot(2,3,4)                                      % 2 EXP 

plot(T4,C4,'ko',T4,(p_pes2(1).*exp(-p_pes2(2).*T4))+(p_pes2(3).*exp(-p_pes2(4).*T4)),'b-')
title('fit con 2 exp WLS')
legend('Data','Fit')
xlabel('tempo [minuti]')
ylabel('concentrazioni [pmol/ml]')


subplot(2,3,5)                                      % 2 EXP SEMILOG
semilogy(T4,C4,'r',T4,(p_pes2(1).*exp(-p_pes2(2).*T4))+(p_pes2(3).*exp(-p_pes2(4).*T4)),'b')
title('semilogaritmica con 2 exp WLS')
legend('Data','Fit')
xlabel('tempo')
ylabel('concentrazioni scala log')


subplot(2,3,6)                                     % 2 EXP COV

scatter(C4,(p_pes2(1).*exp(-p_pes2(2).*T4))+(p_pes2(3).*exp(-p_pes2(4).*T4)),'ko','b')
hold on
plot(tempo,tempo)
axis([0 12 0 12])
title('2 exp WLS')
xlabel('dati')
ylabel('predizioni')


%% ANALISI DEI RESIDUI ESERC 1 e ESERC 3
 
%        2-COMP_WLS                 VS                 2 exp WLS

figure('Name','Analisi dei residui')


% RESIDUI NORMALIZZATI


subplot(2,3,1)                                    % RES NORM 2 COMP
plot(T4,res2comp)
title('residuo norm 2-comp')
hold on
scatter(T4,res2comp,'ko','b')
hold on
yline(0)
yline(-1,'r--')
yline(1,'r--')
yline(-2,'r--')
yline(2,'r--')
legend('Interpolazione','Residui normalizzati')
xlabel('tempo [minuti]')
ylabel('residui norm')


subplot(2,3,4)                                    % RES NORM 2 EXP
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


% RESIDUI su CONCENTRAZIONE


subplot(2,3,2)                                    % RES SU CONC NORM 2 COMP 
scatter(C4,res2comp,'ko','b')
title('RESIDUI IN FUNZIONE DI C norm 2-comp')
hold on
yline(0)
xlabel('concentrazione [pmol/ml]')
ylabel('residui norm')


subplot(2,3,5)                                    % RES SU CONC NORM 2 EXP
scatter(C4,res_pes2,'ko','b')
title('RESIDUI IN FUNZIONE DI C norm 2 exp WLS')
hold on
yline(0)
xlabel('concentrazione [pmol/ml]')
ylabel('residui norm')

% HISTOGRAMMI e GAUSSIANA


subplot(2,3,3)                                       % 2 COMP
histogram(res2comp,10,'Normalization','pdf')
hold on
plot(x_g,Gauss,'r')
title('hist/gauss 2-comp')
legend('hist','gauss')


subplot(2,3,6)                                       % 2 EXP
histogram(res_pes2,10,'Normalization','pdf')
hold on
plot(x_g,Gauss,'r')
title('hist/gauss 2 exp WLS norm')
legend('hist','gauss')
sgtitle('ANALISI RESIDUI 2 COMP CONVERTITO E EXP');



%% MONTECARLO
% CONSIDERO WLS 2 COMP

% valuto la AFFIDABILITA' delle STIME

% OBBIETTIVO:
% GENERAZIONE SINTETICA DI SET DI DATI (STESSA NUMEROSITÀ E ISTANTI DI 
% TEMPO DI QUELLO ORIGINALE) OTTENGO DATI RUMOROSI A PARTIRE DAL VALORE 
% STIMATO, SI RISTIMANO I PARAMETRI POI SI VALUTA LA VARIANZA DELLE STIME

% IL RUMORE SARÀ:
%    NORMALE 
%    MEDIA = 0 
%    VARIANZA = CV*E[X]


giri = 50;     %10000;
pMC = zeros(giri, 4); 

% AD OGNI GIRO GENERO UN SET DI 4 PARAMETRI 

for i=1:giri
    
    % calcolo le misure rumorose, quindi sommo l'errore alle misure

    % poi stimo nuovamente i parametri con i dati affetti da rumore

    % con 'normrnd' che generates a random number from the normal distribution


    % CONSIDERO IL MODELLO LS CON I PARAMETRI COMPARTIMENTALI

    noise = normrnd(0,CV.*C4);
    c_noise = c2comp+noise;

    funMC = @(p_MC) LR(p_MC,T4,c_noise,W_lsqnonlin,dose);
    [p_MC,WRSS_mc,res_mc,exitflag_MC,output_MC,lambda_MC,jacMC]=lsqnonlin(funMC,p2comp);

    pMC(i,:)=p_MC;  % MEMORIZZO IN RIGA I 4 PARAMETRI  [k012,k122,k212,V12]

end


% PRECISIONE delle STIME MONTECARLO
%       CV=SD/E[x]
%       var=SD^2

% se A è una matrice, mean(A) restituisce le medie di ogni colonna di A
E = mean(pMC);
sigMC = std(pMC);
CV_MC = (std(pMC)./E)*100;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% DISTRIBUZIONE dei PARAMETRI MONTECRLO: HISTOGRAMMI e GAUSSIANA

x_gMC = 0:0.001:0.1;
x_gV = 3100:1:4000;

Gauss1 = exp(-(x_gMC-E(1)).^2./(2*sigMC(1)^2))./(sigMC(1)*sqrt(2*pi));
Gauss2 = exp(-(x_gMC-E(2)).^2./(2*sigMC(2)^2))./(sigMC(2)*sqrt(2*pi));
Gauss3 = exp(-(x_gMC-E(3)).^2./(2*sigMC(3)^2))./(sigMC(3)*sqrt(2*pi));
Gauss4 = exp(-(x_gV-E(4)).^2./(2*sigMC(4)^2))./(sigMC(4)*sqrt(2*pi));

figure ('Name','Incertezza sui parametri MONTECARLO')


subplot (2,2,1)
histogram(pMC(:,1),10,'Normalization','pdf')
axis([0.06 0.08 0 300])
title('k01')
hold on
plot(x_gMC,Gauss1,'r')
legend('hist','gauss')

subplot (2,2,2)
histogram(pMC(:,2),10,'Normalization','pdf')
axis([0.05 0.065 0 400])
hold on
plot(x_gMC,Gauss2,'r')
title('k12')
legend('hist','gauss')

subplot (2,2,3)
histogram(pMC(:,3),10,'Normalization','pdf')
axis([0.04 0.08 0 150])
hold on
plot(x_gMC,Gauss3,'r')
title('k21')
legend('hist','gauss')

subplot (2,2,4)
histogram(pMC(:,4),10,'Normalization','pdf')
hold on
plot(x_gV,Gauss4,'r')
title('V')
legend('hist','gauss')

%% ANALISI EDB CON TEST STATISTICI E ANALISI DEI QUARTILI (FACOLTATIVO)
% NORMALITA' - INTERVALLI di CONFIDENZA - PERCENTILI 
% considero MONTECARLO



[h_k01_ks,p_k01_ks] = kstest2(pMC(:,1),Gauss1);  % h=1
[h_k12_ks,p_k12_ks] = kstest2(pMC(:,2),Gauss2);  % h=1
[h_k21_ks,p_k21_ks] = kstest2(pMC(:,3),Gauss3);  % h=1
[h_V_ks,p_V_ks] = kstest2(pMC(:,4),Gauss4);      % h=1


%     #   i dati che ottengo da MC sono normali perche ho la somma
%         con un errore epsilon estratto dalla normale
%     #   le stime invece NON sono una combinazione lineare dei PARAMETRI
%         quindi NON conosco la loro DISTRIBUZIONE
%    MA   so anche che per n-->infinito è asintoticamente normale

%  QUINDI ho due possibiita:
%     1   posso stimare gli intervalli di confidenza dalla normale del 95%
%         sotto hp n>> come E+-2*stdev.
%     2   li posso stimare non conoscendo la distribuzione tramite i
%         PERCENTILI

% ---------------------------------------------------------------------------- %
% 1 -- HP NORMALE n>>

int1 = E+2.*sigMC;
int2 = E-2.*sigMC;

% 2 -- PERCENTILI
perc = prctile(pMC,[5 95]);

% ---------------------------------------------------------------------------- %
% GRAFICO 

%  metodo 1 -- INTERVALLI di CONFIDENZA dei PARAMETRI STIMATI 
figure('Name','Affidabiilta delle stime MONTECARLO')
% figure 6
subplot (4,2,1)
area(x_gMC,Gauss1)
axis([0.06 0.08 0 200])
title('Int confidenza k01')
hold on
xline(E(1),'r--')
xline(int2(1),'r--')
xline(int1(1),'r--')

subplot (4,2,3)
area(x_gMC,Gauss2)
axis([0.05 0.065 0 300])
title('Int confidenza k12')
hold on
xline(E(2),'r--')
xline(int2(2),'r--')
xline(int1(2),'r--')

subplot (4,2,5)
area(x_gMC,Gauss3)
axis([0.04 0.08 0 120])
title('Int confidenza k21')
hold on
xline(E(3),'r--')
xline(int2(3),'r--')
xline(int1(3),'r--')

subplot (4,2,7)
area(x_gV,Gauss4)
title('Int confidenza V')
hold on
xline(E(4),'r--')
xline(int2(4),'r--')
xline(int1(4),'r--')

% ---------------------------------------------------------------------------- %
%  metodo 2 -- PERCENTILI

subplot(4,2,2)
boxplot(pMC(:,1))
title('Percentili k01')

subplot(4,2,4)
boxplot(pMC(:,2))
title('Percentili k12')

subplot(4,2,6)
boxplot(pMC(:,3))
title('Percentili k21')

subplot(4,2,8)
boxplot(pMC(:,4))
title('Percentili V')


%% BOOTSTRAP
% IN GENERE IL METODO DI ESRAZIONE CON REINSERIMENTO TIPICO DEL BOOTSTRAP
% LO SI FA SU MODELLI DI POPOLAZIONE, IN QUESTO CASO INVECE HO UN MODELLO
% CHE DIPENDE DALLA CONCENTRAZIONE, QUINDI NON FACCIO IL REINSERIMENTO SUI
% DATI MA CONSIDRO DI AGGIUNGERE L'ERRORE AI RESIDUI NORMALIZZATI 
% (PERCHÈ SE NON SONO NORMALIZZATI NON POSSO PERMUTARLI)

for i=1:giri
    
    indperm = randi(length(res2comp),30,1); % GENERO INDICI DEL REINSERIMENTO

    resBoo = res2comp(indperm).*(sqrt(sigma)); % ESTRAGGO I RESIDUI NORMALIZZATI

    c_noise_boo = c2comp+resBoo;

    funBOO = @(p_BOO) LR(p_BOO,T4,c_noise_boo,W_lsqnonlin,dose);
    [p_BOO,WRSS_boo,res_boo,exitflag_BOO,output_BOO,lambda_BOO,jacBOO]=lsqnonlin(funBOO,p2comp);

    pBOO(i,:)=p_BOO;

    
end 

% ---------------------------------------------------------------------------- %
% PRECISIONE delle STIME BOOTSTRAP
%       CV=SD/E[x]
%       var=SD^2

% se A è una matrice, mean(A) restituisce le medie di ogni colonna di A

E_BOO = mean(pBOO);
sigBOO = std(pBOO);
CV_BOO = (std(pBOO)./E_BOO)*100;

% PERCENTILI

percBOO = prctile(pBOO,[5 95]);


% ---------------------------------------------------------------------------- %
% DISTRIBUZIONE dei PARAMETRI BOOTSTRAP: HISTOGRAMMI 

figure ('Name','Incertezza sui parametri BOOTSTRAP')
% figure 7

subplot (4,2,1)
histogram(pBOO(:,1),10,'Normalization','pdf')
title('k01')

subplot (4,2,3)
histogram(pBOO(:,2),10,'Normalization','pdf')
title('k12')

subplot (4,2,5)
histogram(pBOO(:,3),10,'Normalization','pdf')
title('k21')

subplot (4,2,7)
histogram(pBOO(:,4),10,'Normalization','pdf')
title('V')

% ---------------------------------------------------------------------------- %
%  PERCENTILI

subplot(4,2,2)
boxplot(pBOO(:,1))
title('Percentili k01')

subplot(4,2,4)
boxplot(pBOO(:,2))
title('Percentili k12')

subplot(4,2,6)
boxplot(pBOO(:,3))
title('Percentili k21')

subplot(4,2,8)
boxplot(pBOO(:,4))
title('Percentili V')


%% REPORT CON DATI SIGNIFICATIVI

fid = fopen('eserc4modelli.txt','w');

fprintf(fid,'\n DATI RILEVANTI ESERCITAZIONE 4 \n');

fprintf(fid,'\nmodello      A        alpha          B        beta');
fprintf(fid,'\n%s\t%f\t%f\t%f\t%f\t','2-comp',A_2comp,Alpha_2comp,B_2comp,Beta_2comp);
fprintf(fid,'\n%s\t%f\t%f\t%f\t%f\t','2 exp',A_stim,Alpha_stim,B_stim,Beta_stim);

fprintf(fid,'\n\n\n\n CONFRONTO PARAMETRI DATI e STIMATI \n');
fprintf(fid,'\n#        K01           K12         K21         V');
fprintf(fid,'\n%s\t%f\t%f\t%f\t%f\t','dati',p2comp(1),p2comp(2),p2comp(3),p2comp(4));
fprintf(fid,'\n%s\t\t%f\t%f\t%f\t%f\t','2-c',p2_comp(1),p2_comp(2),p2_comp(3),p2_comp(4));
fprintf(fid,'\n%s\t\t%f\t%f\t%f\t%f\t','MC',E(1),E(2),E(3),E(4));
fprintf(fid,'\n%s\t\t%f\t%f\t%f\t%f\t','BOO',E_BOO(1),E_BOO(2),E_BOO(3),E_BOO(4));



fprintf(fid,'\n\n\n\n CV delle STIME Montecarlo, Bootstrap e WLS-2comp \n');
fprintf(fid,'\n#      K01         K12         K21         V');
fprintf(fid,'\n%s\t%f\t%f\t%f\t%f\t','MC',CV_MC(1),CV_MC(2),CV_MC(3),CV_MC(4));
fprintf(fid,'\n%s\t%f\t%f\t%f\t%f\t','BOO',CV_BOO(1),CV_BOO(2),CV_BOO(3),CV_BOO(4));
fprintf(fid,'\n%s\t%f\t%f\t%f\t%f\t','2-c',CV_2comp(1),CV_2comp(2),CV_2comp(3),CV_2comp(4));

fprintf(fid,'\n\n\n\n INTERVALLI di CONFIDENZA stime Montecarlo \n');
fprintf(fid,'\n#      SX           DX');
fprintf(fid,'\n%s\t%f\t%f\t','K01',int2(1),int1(1));
fprintf(fid,'\n%s\t%f\t%f\t','K01',int2(2),int1(2));
fprintf(fid,'\n%s\t%f\t%f\t','K01',int2(3),int1(3));
fprintf(fid,'\n%s\t%f\t%f\t','V',int2(4),int1(4));

fprintf(fid,'\n\n\n\n PERCENTILI stime Montecarlo e Bootstrap \n');
fprintf(fid,'\n#        5 MC            95 MC          5 BOO           95 BOO');
fprintf(fid,'\n%s\t\t%f\t\t%f\t\t%f\t\t%f\t','K01',perc(1,1),perc(2,1),percBOO(1,1),percBOO(2,1));
fprintf(fid,'\n%s\t\t%f\t\t%f\t\t%f\t\t%f\t','K12',perc(1,2),perc(2,2),percBOO(1,2),percBOO(2,2));
fprintf(fid,'\n%s\t\t%f\t\t%f\t\t%f\t\t%f\t','K21',perc(1,3),perc(2,3),percBOO(1,3),percBOO(2,3));
fprintf(fid,'\n%s\t\t%f\t\t%f\t\t%f\t\t%f\t','V',perc(1,4),perc(2,4),percBOO(1,4),percBOO(2,4));

fprintf(fid,'\n\n\n\n\n');
fclose(fid);


