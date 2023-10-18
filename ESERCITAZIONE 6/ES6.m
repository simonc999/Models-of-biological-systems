%ESERCITAZIONE 6 - REGRESSIONE LINEARE

%OBIETTIVO: STIMARE I PARAMETRI DEL MODELLO CHE DESCRIVE LA CINETICA DEL
%C-PEPTIDE (Y) A PARTIRE DA REGRESSORI ANTROPOMETRICI CHE CARATTERIZZANO IL
%SOGGETTO, COME SUGGERITO DALL'APPROCCIO DI VAN CAUTER. VAN CAUTER PROPOSE
%DI DEFINIRE UN MODELLO CHE SFRUTTASSE LA LIMITATA VARIABILITÀ DELLA
%CINETICA NEI DIVERSI SOGGETTI, NOTI ALCUNI PARAMETRI ANTROPOMETRICI. I
%PARAMETRI CHE SI VOGLIONO STIMARE NON SONO DIRETTAMENTE QUELLI DEL MODELLO
%COMPARTIMENTALE, MA PARAMETRI (BETA) DERIVATI CHE MIGLIORANO LA RELAZIONE
%LINEARE FRA COVARIATE E PARAMETRI.

% PARALLELISMO CON MODELLI DI POPOLAZIONE: L'EFFETTO RANDOM (ETA) È QUELLO
% CHE NELLA REGRESSIONE ABBIAMO CHIAMATO EPSILON (ERRORE RESIDUO). I TETA
% DEI MODELLI DI POPOLAZIONE SONO QUELLI CHE QUI ABBIAMO CHIAMATO BETA. 

% IN QUESTA ESERCITAZIONE NON ABBIAMO STIMATO I PARAMETRI DI POPOLAZIONE MA
% ABBIAMO STIMATO I VALORI A PARTIRE DAI PARAMETRI INDIVIDUALI CHE A SUA VOLTA 
% QUALCUNO AVEVA STIMATO A PARTIRE DALLE MISURE DI C PEPTIDE. 

clear 
close all
clc
run dati_vc

Dati=table(dati(:,1),dati(:,2),dati(:,3),dati(:,4),dati(:,5),dati(:,6),dati(:,7),dati(:,8),dati(:,9),dati(:,10),dati(:,11));
Dati.Properties.VariableNames={'Statodisalute','Sesso','Eta','Peso','Altezza','BMI','BSA','Vd','Emivitacorto','Emivitalungo','Frazione'};
grandezze=[Dati.Vd,Dati.Emivitacorto,Dati.Emivitalungo,Dati.Frazione];

%% 1) %%%%%%%%------------  ANALISI MONOVARIATA  -------------%%%%%%%%%%%%%

% grafico i dati raccolti in fuzione delle singole covariate.

% PER LE VARIABILI CATEGORICHE POSSO FARE BOXPLOT, OGNI VALORE HO UNA 
% DISTRIBUZIONE DEI VALORI VD, EMIVITA CORTO E LUNGO E FRACTION: 
% MEDIANA, RANGE CHE COMPRENDE DAL 25 AL 75% PERCENTILE E EVENTUALI 
% OUTLIERS. 
% HO DUE FIGURE CON I BOXPLOT PER SESSO E SALUTE
% PER LE VARIABILI CONTINUE POSSO RAPPRESENTARE I DATI IN FUNZIONE DELLA
% VARIABILE, UNA FIG PER OGNI VARIABILE CONTINUA. SI EVIDENZIA CHE UNA
% VARIABILE NON È SUFFICIENTE A SPIEGARE LA VARIABILITÀ DEI VALORI SECONDO
% RELAZIONE LIN: PUNTI NON ALLINEATI LUNGO RETTA, DISPOSTI CASUALMENTE.


for i=1:7
    figure(i);
    for j=1:4
        if i<3 % HA SENSO SOLO PER STATO E SESSO QUESTO BOXPLOT VARIABILI 
               % CATEGORICHE

            subplot(2,2,j);
            boxplot(grandezze(:,j),dati(:,i)); 
            % I PIÙ SONO GLI OULIER FUORI
            % DA 25 E 75 PERCENTILE CHE SONO 1/4 DEI DATI E 3/4 DEI DATI

        if i==1
            xlabel({'0=normale, 1=obeso, 2=diabetico'});
        end
        if i==2
            xlabel({'0=M, 1=F'});
        end
        ylabel(Dati.Properties.VariableNames{j+7});
        sgtitle(Dati.Properties.VariableNames{i});
    else
        subplot(2,2,j);
        scatter(dati(:,i),grandezze(:,j),'*');
        xlabel(Dati.Properties.VariableNames{i});
        ylabel(Dati.Properties.VariableNames{j+7});
        sgtitle(Dati.Properties.VariableNames{i});
        end
    end
end

% VARIABILI DUMMY, HO VARIABILI CONTINUE, 1 O 0 PER MODELLI CON VARIABILI 
% CATEGORICHE
%            0           1         2

% CODIFICA: 00 NORMALE, 01 OBESO, 10 DIABETICO
% DIVISIONE IN DUE REGRESSORI DELLE VARIABILI DUMMY PER LO STATO DI SALUTE

% CONVERTO LA VARIABILE STATO DI SALUTE IN DUMMY VARIABLE

stato_salute=dati(:,1);
for j=1:length(stato_salute)
    if stato_salute(j)==0
       X1(j)=0;
       X2(j)=0;
    end
    if stato_salute(j)==1
      X1(j)=0;
      X2(j)=1;
    end
    if stato_salute(j)==2
       X1(j)=1;
       X2(j)=0;      
    end
end


%% 2) %%%%%%%%-------TROVO MODELLO OTTIMALE REGRESSIONE-------%%%%%%%%%%%%%


% successivamente sulla base dei dati disponibili, 
% 1 trovare il miglior modello di regressione, 
% 2 determinare i coefficenti della regressione di  tale modello 
% 3 con i loro cv e 4 gli intervalli di confidenza.
% 5 verificare la significativit`a statistica dei regressori inclusi nel 
% modello.

% PER LE VARIABILI CATEGORICHE NON DICOTOMICHE USO DUMMY VARIABLES, CON N
% CATEGORIE USO N-1 DUMMY: 00 NORMALE, 01 OBESO, 10 DIABETICO.
% SE FOSSERO STATE 4 000 001 010 100
% DEFINISCO A MATRICE REGRESSORI CONTENENTE I SINGOLI E LE COMBINAZIONI DI
% VARIABILE CATEGORICA E CONTINUA, DI VARIABILE CATEGORICA E CATEGORICA. 
% PER OGNI GRANDEZZA DA STIMARE RIFACCIO A E RIPETO ALGORITMO.

% PER DETERMINARE I REGRESSORI DA USARE VIENE UTILIZZATO FORWARD STEPWISE.
% MODELLO DI REGRESSIONE LINEARE MULTIPLA PER VOLUME DI DISTRIBUZIONE 


%-BSA e età  per il Vd
 %-BMI e età per emivita corto
 %-eta e altezza per emivita lungo
% -BMI per fraction

% CELLARRAY CONTENTE I REGRESSORI. 18 CELLE: NELLA PRIMA CELLA SI TROVANO
% LE DUE VARIABILI DUMMY. LE COMBINAZIONI PREVEDONO IL PRODOTTO TRA 
% VARIABILI E QUELLE CONTINUE

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


A={[X1', X2'],dati(:,2),dati(:,3),dati(:,4),dati(:,5),dati(:,6),dati(:,7),... %% 1--
dati(:,2).*[X1' X2'] ,dati(:,2).*dati(:,3), dati(:,2).*dati(:,4),...          %% 2--
    dati(:,2).*dati(:,5),dati(:,2).*dati(:,6), dati(:,2).*dati(:,7),...       %% 2--
[X1' X2'].*dati(:,3), [X1' X2'].*dati(:,4), [X1' X2'].*dati(:,5),...          %% 3--
    [X1' X2'].*dati(:,6), [X1' X2'].*dati(:,7)};                              %% 3--

% 1--ANALIZZO PRIMA I REGRESSORI SINGOLI INSERENDO NELLO STATO DI SALUTE LE
% DUMMY VARIABLES  ---> 7 REGRESSORI

% ...
    
% 2--COMBINO PER LA VARIABILE CATEGORICA "SESSO" LE ALTRE VARIABILI
%                  ---> 6 REGRESSORI
% ...

    
% 3--COMBINO PER LA VARIABILE CATEGORICA "STATO DI SALUTE" LE ALTRE VARIABILI
%                  ---> 5 REGRESSORI
  
% A = MATRICE DEI REGRESSORI 
% NON FISSO LA PENDENZA E LA INTERCETTA DELLA RETTA IN QUESTO MODO.

%% MODELLO DI REGRESSIONE LINEARE MULTIPLA PER VOLUME DI DISTRIBUZIONE

Y=dati(:,8);
n=length(Y);

X0=[ones(n,1)]; % CONSIDERO BETA0

tot_reg=size(A,2); % RICAVO IL NUMERO DI COLONNE DELLA MATRICE DELLA 
                   % COMBINAZIONE DEI REGRESSORI

reg=tot_reg; % INIZIALMENTE DEVONO ESSERE PROVATI TUTTI, UNO PER VOLTA. 
             % IL VALORE REG DIMINUISCE DURANTE L'ALGORITMO, QUANDO 
             % INSERISCO UN NUOVO REGRESSORE. IN TOTALE SONO 18.

             
R2adj_rif=0; % INIZIALMENTE FISSO R2 ADJ DI CONFRONTO PARI A ZERO

for i=1:tot_reg  
    R2adj=[];  

    for j=1:reg
        X=[X0 A{:,j}]; % AD OGNI ITERAZIONE DEL CICLO INTERNO X VIENE 
                       % AGGIORNATA CON UNA COLONNA ALLA VOLTA DELLA 
                       % MATRICE DEI REGRESSORI IN MODO CHE OGNI REGRESSORE
                       % SIA VALUTATO.
                       % 
                       % LA MATRICE X HA TANTE RIGHE QUANTI SONO I SOGGETTI
                       % E TANTE COLONNE QUANTI SONO I REGRESSORI
                       % CONSIDERATI + 1 RIFERITO ALLA COLONNA IDENTITA'
                       % RELATIVA AI BETA0

        betaS=(inv(X'*X))*X'*Y; % FORMULA STIME BETA

        yP=X*betaS; % CALCOLO PREDIZIONI (VETTORE DELLE PREDIZIONI USANDO 
                    % IL VETTORE DEI BETA STIMATI)
                    % STIMIAMO IL VOLUME SENZA ERRORE PER POI PARAGONARLO
                    % ALLE MISURE E TROVARE L'ERRORE STESSO

        SST=sum((Y-mean(Y)).^2); % somma dei quadrati totale (VAR TOTALE)
        SSE=sum((Y-yP).^2);      % somma dei quadrati dell'errore (NON SPIEGATA)

        % RICAVO L'INDICE R2 CHE E' IL RAPPORTO TRA DEVIANZA SPIEGATA /
        % TOTALE

        R2=(SST-SSE)/SST;

        k=length(betaS)-1; % N REGRESSORI USATI, CALCOLANDO K POSSO TROVARE 
                           % R2 ADJ


        R2adj(j)=R2-(k/(n-k-1))*(1-R2); 

                           % SALVO GLI R2 ADJ DI OGNI REGRESSORE PER 
                           % SCEGLIERE IL MAGGIORE
    end

    % AD OGNI CICLO ESTERNO SALVO AD OGNI CICLO IL MASSIMO R2ADJ E IL 
    % SUO INDICE

    [R2adj_max(i), idxMax(i)]=max(R2adj); 

    if ne(i,1) % AL PRIMO GIRO IL TERMINE DI CONFRONTO RIMANE 0
        R2adj_rif=R2adj_max(i-1);
    end
    if R2adj_max(i) <= R2adj_rif % SE IL MASSIMO CORRENTE È MINORE DEL 
                                 % MASSIMO PRECEDENTE MI FERMO
        X=X0;
        break;
    else
        X0=[X0 A{:,idxMax(i)}]; % ALTRIMENTI AGGIORNO LA MATRICE X0 CON IL 
                                % REGRESSORE CORRISPONDENTE A R2 ADJ 
                                % MASSIMO E RIMUOVO LA COLONNA DEL
                                % REGRESSORE SCELTO DALLA MATRICE A
        A{:,idxMax(i)}=[];
        reg=reg-1;
    end
end
    
% R2ADJ È IL VETTORE DELL'R^2 ADJ DEL VOLUME DI DISTRIBUZIONE

% R2ADJ MAX È IL VETTORE CHE CONTIENE I VALORI MASSIMI DELL'R^2 ADJ 

% CON LA MATRICE X CHE SPIEGA IL MASSIMO DELLA DEVIANZA CALCOLO I BETA, LE
% STIME E I RESIDUI

betaSF=(inv(X'*X))*X'*Y; % BETA NON DIPENDE DA ALFA (COSTANTE MOLTIPLICATIVA DOVUTA ALLA SIGMA)

yPF=X*betaSF;

residui=Y-yPF; % CALCOLO I RESIDUI

figure(8);

subplot(1,2,1);                                                 % RESIDUI
plot(yPF,residui,'*');
xlabel('Vd');
ylabel('residui');
title('Residui');

subplot(1,2,2);                                                 % QQ RES
% (*QQ-verifico-distr-normale*)
qqplot(residui); % VERIFICO SE SEGUONO UNA DISTRIBUZIONE NORMALE

sgtitle('Analisi residui Vd');

SSE=sum((Y-yPF).^2);
p=length(betaSF)-1;  % NUMERO DEI REGRESSORI
R2FIN=(SST-SSE)/SST; % R2 FINALE DEL MODELLO 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%% PRECISIONE DEI PARAMETRI STIMATI: DEVO UTILIZZARE LA STIMA DELLA VARIANZA
% DELL'ERRORE: NON CONOSCO SIGMA_V DEVO STIMARLA, SIGMA_V=MATRICE DIAGONALE
% CON DIAGONALE S^2(VARIANZA).

s2=SSE/(n-p-1);  % VARIANZA NON SPIEGATA


% EQUIVALENTE DELLA SIGMA V
sigmavS=s2.*eye(n);   % L'ASSUNZIONE SULL'ERRORE PREVEDE CHE LA MATRICE 
                      % SIGMAV SIA DIAGONALE E CHE LA VARIANZA SIA 
                      % COSTANTE, USO LA VARIANZA STIMATA 


% EQUIVALENTE DELLA SIGMA P
sigma_betaSF=inv(X'*inv(sigmavS)*X); % MATRICE DELLE VARIANZE DELLE STIME, 
                                     % FORMULA CHIUSA

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% CALCOLO DEI CV E DEGLI INTERVALLI DI CONFIDENZA DEI PARAMETRI BETA


CV=zeros(3,1);

alpha=0.05;

t=tinv(1-alpha/2,n-p-1);

intervalli=cell(3,1);

for i=1:length(betaSF)

    CV(i)=sqrt(sigma_betaSF(i,i))/abs(betaSF(i)); 
    intervalli{i}=[betaSF(i)-t*sqrt(sigma_betaSF(i,i)) betaSF(i)+t*sqrt(sigma_betaSF(i,i))];

end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% VALUTAZIONE DELLA SIGNIFICATIVITÀ STATISTICA DEI REGRESSORI: USO TEST F
% TRA MODELLO COMPLETO E RIDOTTO, CON BACKWARS STEPWISE.

% 1) VERIFICO SIGNIFICATIVITÀ DEL REGRESSORE ETA'

% H0: L'ETA' NON È REGRESSORE SIGNIFICATIVO PER VOLUME DI DISTRIBUZIONE, 
% ALPHA=0.05

% RIDUCO IL MODELLO

X_red=[X(:,1),X(:,2)]; % TOLGO LA COLONNA RELATIVA ALL'ULTIMO REGRESSORE

betaSF_red=(inv(X_red'*X_red))*X_red'*Y;

yP_red=X_red*betaSF_red;

SSred=sum((Y-yP_red).^2);

k=1; % CONSIDERO UN SOLO REGRESSORE

F=((SSred-SSE)/(p-k))/(SSE/(n-p-1)); 

pvalue=fcdf(F,p-k,n-p-1,'upper');   % pvalue < alpha ALLORA RIFIUTO H0

% 2) VERIFICO SIGNIFICATIVITÀ DEL REGRESSORE BSA
% H0: BSA NON È REGRESSORE SIGNIFICATIVO PER VOLUME DI DISTRIBUZIONE

% RIDUCO IL MODELLO

X_red2=[X(:,1),X(:,3)]; % QUESTA VOLTA TOLGO LA COLONNA RELATIVA AL PRIMO 
                        % REGRESSORE

betaSF_red2=(inv(X_red2'*X_red2))*X_red2'*Y;
yP_red2=X_red2*betaSF_red2;
SSred2=sum((Y-yP_red2).^2);
k=1; % CONSIDERO UN SOLO REGRESSORE

F2=((SSred2-SSE)/(p-k))/(SSE/(n-p-1)); 

pvalue2=fcdf(F2,p-k,n-p-1,'upper');

% H0 VIENE RIFIUTATA

%% MODELLO DI REGRESSIONE LINEARE MULTIPLA PER EMIVITA CORTO

A={[X1' X2'],dati(:,2),dati(:,3),dati(:,4),dati(:,5),dati(:,6),dati(:,7),...
    dati(:,2).*[X1' X2'],dati(:,2).*dati(:,3),dati(:,2).*dati(:,4),dati(:,2).*dati(:,5),dati(:,2).*dati(:,6), dati(:,2).*dati(:,7),...
    [X1' X2'].*dati(:,3),[X1' X2'].*dati(:,4),[X1' X2'].*dati(:,5),[X1' X2'].*dati(:,6),[X1' X2'].*dati(:,7)};

YE=dati(:,9);
n=length(YE);
X0=ones(n,1);

tot_reg=size(A,2);
reg=tot_reg;
R2adj_rif=0;

for i=1:tot_reg

    R2adjE=[];
    for j=1:reg
        XE=[X0 A{:,j}]; 
        betaS=(inv(XE'*XE))*XE'*YE;
        yP=XE*betaS;
        SST=sum((YE-mean(YE)).^2);
        SSE=sum((YE-yP).^2);
        R2=(SST-SSE)/SST;
        k=length(betaS)-1; 
        R2adjE(j)=R2-(k/(n-k-1))*(1-R2); 
    end
    [R2adj_maxE(i) idxMaxE(i)]=max(R2adjE); 

    if ne(i,1) 
        R2adj_rif=R2adj_maxE(i-1);
    end

    if R2adj_maxE(i) <= R2adj_rif 
        XE=X0;
        break;
    else
        X0=[X0 A{:,idxMaxE(i)}];
        A{:,idxMaxE(i)}=[];
        reg=reg-1;
    end
end

% IL MODELLO DI REGRESSIONE TROVATO PREVEDE REGRESSORI BMI E ETA

betaSFE=(inv(XE'*XE))*XE'*YE;
yPF=XE*betaSFE;
residuiE=YE-yPF;


figure(9);

subplot(1,2,1);                                                 % RESIDUI
plot(yPF,residuiE,'*');
xlabel('Emivita Corto');
ylabel('residui');
title('Residui');

subplot(1,2,2);                                                 % QQ RES
qqplot(residuiE);
sgtitle('Analisi dei residui per Emivita Corto');

SSE=sum((YE-yPF).^2);
p=length(betaSFE)-1; 
R2EFIN=(SST-SSE)/SST;

% PRECISIONE

s2E=SSE/(n-p-1);
sigmavSE=s2E.*eye(n); 
sigma_betaSFE=inv(XE'*inv(sigmavSE)*XE); 

% CALCOLO DEI CV E DEGLI INTERVALLI DI CONFIDENZAù
CVE=zeros(3,1);
alpha=0.05;
t=tinv(1-alpha/2,n-p-1);
intervalliE=cell(3,1);

for i=1:length(betaSFE)
    CVE(i)=sqrt(sigma_betaSFE(i,i))/abs(betaSFE(i)); 
    intervalliE{i}=[betaSFE(i)-t*sqrt(sigma_betaSFE(i,i)) betaSFE(i)+t*sqrt(sigma_betaSFE(i,i))];
end

% VALUTAZIONE DELLA SIGNIFICATIVITÀ STATISTICA DEI REGRESSORI
% 1) VERIFICO SIGNIFICATIVITÀ DEL REGRESSORE ETA
% H0: L'ETA NON È REGRESSORE SIGNIFICATIVO, ALPHA=0.05

% RIDUCO IL MODELLO

X_red=[XE(:,1),XE(:,2)]; % TOLGO LA COLONNA RELATIVA AL REGRESSORE ETA
betaSF_red=(inv(X_red'*X_red))*X_red'*YE;
yP_red=X_red*betaSF_red;
SSred=sum((YE-yP_red).^2);
k=1; % CONSIDERO UN SOLO REGRESSORE
FE=((SSred-SSE)/(p-k))/(SSE/(n-p-1)); 
pvalueE=fcdf(FE,p-k,n-p-1,'upper');

% H0 NON VIENE RIFIUTATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% 2) VERIFICO SIGNIFICATIVITÀ DEL REGRESSORE BMI
% H0: BMI NON È REGRESSORE SIGNIFICATIVO

 
%RIDUCO IL MODELLO
X_red2=[XE(:,1),XE(:,3)]; % TOLGO LA COLONNA RELATIVA AL REGRESSORE BMI
betaSF_red2=(inv(X_red2'*X_red2))*X_red2'*YE;
yP_red2=X_red2*betaSF_red2;
SSred2=sum((YE-yP_red2).^2);
k=1; % CONSIDERO UN SOLO REGRESSORE
FE2=((SSred2-SSE)/(p-k))/(SSE/(n-p-1)); 
pvalueE2=fcdf(FE2,p-k,n-p-1,'upper');

% H0 RIFIUTATA


%% MODELLO DI REGRESSIONE LINEARE MULTIPLA PER EMIVITA LUNGO


A={[X1' X2'],dati(:,2),dati(:,3),dati(:,4),dati(:,5),dati(:,6),dati(:,7),...
    dati(:,2).*[X1' X2'],dati(:,2).*dati(:,3),dati(:,2).*dati(:,4),dati(:,2).*dati(:,5),dati(:,2).*dati(:,6), dati(:,2).*dati(:,7),...
    [X1' X2'].*dati(:,3),[X1' X2'].*dati(:,4),[X1' X2'].*dati(:,5),[X1' X2'].*dati(:,6),[X1' X2'].*dati(:,7)};

tot_reg=size(A,2);
YEl=dati(:,10);
n=length(YEl);

X0=[ones(n,1)];
reg=tot_reg; 
R2adj_rif=0;

for i=1:tot_reg

    R2adjEl=[];
    for j=1:reg
        XEl=[X0 A{:,j}]; 
        betaS=(inv(XEl'*XEl))*XEl'*YEl;
        yP=XEl*betaS;
        SST=sum((YEl-mean(YEl)).^2);
        SSE=sum((YEl-yP).^2);
        R2=(SST-SSE)/SST;
        k=length(betaS)-1;
        R2adjEl(j)=R2-(k/(n-k-1))*(1-R2); 
    end
    [R2adj_maxEl(i) idxMaxEl(i)]=max(R2adjEl);

    if ne(i,1) 
        R2adj_rif=R2adj_maxEl(i-1);
    end

    if R2adj_maxEl(i) <= R2adj_rif 
        XEl=X0;
        break;
    else

        X0=[X0 A{:,idxMaxEl(i)}];
        A{:,idxMaxEl(i)}=[];
        reg=reg-1;
    end
end

% IL MODELLO DI REGRESSIONE TROVATO PREVEDE REGRESSORI ETA E ALTEZZA

betaSFEl=(inv(XEl'*XEl))*XEl'*YEl;
yPF=XEl*betaSFEl;
residuiEl=YEl-yPF;

figure(10);

subplot(1,2,1);                                                 % RESIDUI
plot(yPF,residuiEl,'*');
xlabel('Emivita Lungo');
ylabel('residui');
title('Residui');

subplot(1,2,2);                                                 % QQ RES
qqplot(residuiEl);
sgtitle('Analisi dei residui per Emivita Lungo');

SSE=sum((YEl-yPF).^2);
p=length(betaSFEl)-1; 
R2ElFIN=(SST-SSE)/SST;

% PRECISIONE DEI PARAMETRI

s2El=SSE/(n-p-1);
sigmavSEl=s2El.*eye(n);
sigma_betaSFEl=inv(XEl'*inv(sigmavSEl)*XEl);


% CV E INTERVALLI DI CONFIDENZA

CVEl=zeros(3,1);
alpha=0.05;
t=tinv(1-alpha/2,n-p-1);
intervalliEl=cell(3,1);
for i=1:length(betaSFEl)
    CVEl(i)=sqrt(sigma_betaSFEl(i,i))/abs(betaSFEl(i)); 
    intervalliEl{i}=[betaSFEl(i)-t*sqrt(sigma_betaSFEl(i,i)) betaSFEl(i)+t*sqrt(sigma_betaSFEl(i,i))];
end

% SIGNIFICATIVITA' STATISTICA REGRESSORI

% 1) VERIFICO SIGNIFICATIVITÀ DEL REGRESSORE ALTEZZA
% H0: L'ALTEZZA NON È REGRESSORE SIGNIFICATIVO, ALPHA=0.05

% RIDUCO IL MODELLO
X_red=[XEl(:,1),XEl(:,2)]; % TOLGO LA COLONNA RELATIVA AL REGRESSORE ALTEZZA
betaSF_red=(inv(X_red'*X_red))*X_red'*YEl;
yP_red=X_red*betaSF_red;
SSred=sum((YEl-yP_red).^2);
k=1;% CONSIDERO UN SOLO REGRESSORE

FEl=((SSred-SSE)/(p-k))/(SSE/(n-p-1)); 
pvalueEl=fcdf(FEl,p-k,n-p-1,'upper');

% H0 RIFIUTATA

% 2) VERIFICO SIGNIFICATIVITÀ DEL REGRESSORE ETA
% H0: ETA NON È REGRESSORE SIGNIFICATIVO

% RIDUCO IL MODELLO

X_red2=[XEl(:,1),XEl(:,3)]; % TOLGO LA COLONNA RELATIVA AL REGRESSORE ETA
betaSF_red2=(inv(X_red2'*X_red2))*X_red2'*YEl;
yP_red2=X_red2*betaSF_red2;
SSred2=sum((YEl-yP_red2).^2);
k=1;% CONSIDERO UN SOLO REGRESSORE

FEl2=((SSred2-SSE)/(p-k))/(SSE/(n-p-1)); 
pvalueEl2=fcdf(FEl2,p-k,n-p-1,'upper');

% H0 RIFIUTATA

%% MODELLO DI REGRESSIONE LINEARE MULTIPLA PER FRACTION 
% (RAPPORTO TRA I DUE ESPONENZIALI DELLA RISPOSTA IMPULSIVA A1/(A1+A2) 
% DOVE A1 E A2 SONO I COEFFICIENTI DEGLI ESPONENZIALI RELATIVI 
% RISPETTIVAMENTE ALLE COSTANTI DI TEMPO CORTA E LUNGA)


A={[X1' X2'],dati(:,2),dati(:,3),dati(:,4),dati(:,5),dati(:,6),dati(:,7),...
    dati(:,2).*[X1' X2'],dati(:,2).*dati(:,3),dati(:,2).*dati(:,4),dati(:,2).*dati(:,5),dati(:,2).*dati(:,6), dati(:,2).*dati(:,7),...
    [X1' X2'].*dati(:,3),[X1' X2'].*dati(:,4),[X1' X2'].*dati(:,5),[X1' X2'].*dati(:,6),[X1' X2'].*dati(:,7)};

tot_reg=size(A,2);
YF=dati(:,11);
n=length(YF);

X0=[ones(n,1)];
reg=tot_reg; 
R2adj_rif=0;

for i=1:tot_reg

    R2adjF=[];
    for j=1:reg
        XF=[X0 A{:,j}]; 
        betaS=(inv(XF'*XF))*XF'*YF;
        yP=XF*betaS;
        SST=sum((YF-mean(YF)).^2);
        SSE=sum((YF-yP).^2);
        R2=(SST-SSE)/SST;
        k=length(betaS)-1;
        R2adjF(j)=R2-(k/(n-k-1))*(1-R2); 
    end
    [R2adj_maxF(i) idxMaxF(i)]=max(R2adjF); 

    if ne(i,1)
        R2adj_rif=R2adj_maxF(i-1);
    end

    if R2adj_maxF(i) <= R2adj_rif 

        XF=X0;
        break;
    else

        X0=[X0 A{:,idxMaxF(i)}];
        A{:,idxMaxF(i)}=[];
        reg=reg-1;

    end
end

% IL MODELLO DI REGRESSIONE TROVATO PREVEDE REGRESSORI BMI 


betaSFF=(inv(XF'*XF))*XF'*YF;
yPF=XF*betaSFF;
residuiF=YF-yPF;

figure(11);

subplot(1,2,1);                                                 % RESIDUI
plot(yPF,residuiF,'*');
xlabel('Fraction');
ylabel('residui');
title('Residui');

subplot(1,2,2);
qqplot(residuiF);                                               % QQ RES
sgtitle('Analisi dei residui per Fraction');

SSE=sum((YF-yPF).^2);
p=length(betaSFF)-1; 
R2FFIN=(SST-SSE)/SST;

% PRECISIONE PARAMETRI STIMATI

s2F=SSE/(n-p-1);
sigmavSF=s2F.*eye(n);
sigma_betaSFF=inv(XF'*inv(sigmavSF)*XF); 

% CALCOLO CV E INT DI CONFIDENZA

CVF=zeros(2,1);
alpha=0.05;
t=tinv(1-alpha/2,n-p-1);
intervalliF=cell(3,1);
for i=1:length(betaSFF)
    CVF(i)=sqrt(sigma_betaSFF(i,i))/abs(betaSFF(i)); 
    intervalliF{i}=[betaSFF(i)-t*sqrt(sigma_betaSFF(i,i)) betaSFF(i)+t*sqrt(sigma_betaSFF(i,i))];
end

% VALUTAZIONE DELLA SIGNIFICATIVITÀ STATISTICA DEI REGRESSORI
% 1) VERIFICO SIGNIFICATIVITÀ DEL REGRESSORE BMI
% H0: BMI NON È REGRESSORE SIGNIFICATIVO, ALPHA=0.05

% RIDUCO IL MODELLO

X_red=[XF(:,1)]; % TOLGO LA COLONNA RELATIVA AL REGRESSORE BMI
betaSF_red=(inv(X_red'*X_red))*X_red'*YF;
yP_red=X_red*betaSF_red;
SSred=sum((YF-yP_red).^2);
k=0; % CONSIDERO 0 REGRESSORI
FF=((SSred-SSE)/(p-k))/(SSE/(n-p-1)); 
pvalueF=fcdf(FF,p-k,n-p-1,'upper');

% H0
% RIFIUTATA




%% 3)Fare la predizione di vostri parametri della cinetica del CP usando i
% vostri dati (sesso, peso, altezza,...). In particolare, fate la 
% predizione puntuale e determinate intervalli di confidenza dei
% parametri della cinetica (nelle 3 possibili parametrizzazioni).

% USO I BETA TROVATI, PER EMIVITA CORTO ELIMINA IL SUO BETA

% PRIMA PARAMETRIZZAZIONE: VD, EMIVITA CORTO, EMIVITA LUNGO, FRACTION 

%% PREDIZIONE DEL VD: REGRESSORI BSA E ETA


alpha=0.05;
n=length(Dati.Vd);
t=tinv(1-alpha/2,n-3);

eta=21;                                    %%%%%%%%% ! %%%%%%%%%%
altezza=1.83;
peso=80;

BSA=0.20247*(altezza^0.725)*(peso^0.425);  %%%%%%%%% ! %%%%%%%%%% 
% RICAVATO CON FORMULA

Xi=[1 BSA eta]; % SCRIVO LA MATRICE DEI REGRESSORI

% Y = X*BETA

Vd=Xi*betaSF; % RICAVO Y (VOLUME DI DISTRIBUZIONE)

s2Vd=Xi*(inv(X'*inv(sigmavS)*X))*Xi'; % VARIANZA DELLE PREDIZIONI STIMATE

intervalloVd=[Vd-t*sqrt(s2Vd) Vd+t*sqrt(s2Vd)]; % TROVO INTERVALLO DI Y^


%% PREDIZIONE EMIVITA CORTO: REGRESSORI BMI E ETA

BMI=peso/(altezza^2);                                %%%%%%%%% ! %%%%%%%%%%
XiE=[1 BMI eta];
emivita_corto=XiE*betaSFE;
s2E=XiE*(inv(XE'*inv(sigmavSE)*XE))*XiE';
intervalloE=[emivita_corto-t*sqrt(s2E) emivita_corto+t*sqrt(s2E)];

%% PREDIZIONE EMIVITA LUNGO: REGRESSORI ETA E ALTEZZA

XiEl=[1 eta altezza];
emivita_lungo=XiEl*betaSFEl;
s2El=XiEl*(inv(XEl'*inv(sigmavSEl)*XEl))*XiEl';
intervalloEl=[emivita_lungo-t*sqrt(s2El) emivita_lungo+t*sqrt(s2El)];

%% PREDIZIONE FRACTION

XiF=[1 BMI];
fraction=XiF*betaSFF;
s2F=XiF*(inv(XF'*inv(sigmavSF)*XF))*XiF';
intervalloF=[fraction-t*sqrt(s2F) fraction+t*sqrt(s2F)];



% DISTRIBUZIONE: SO CHE VD, EMIVITA CORTO, EMIVITA LUNGO E FRACTION SONO
% DISTRIBUITI COME T STUDENT; DALLA STATISTICA T STANDARDIZZATA, DEVO
% TORNARE AL VALORE AGGIUNGENDO LA MEDIA E MOLTIPLICANDO PER LA DEVIAZIONE
% STANDARD


n=207;
dim=[100 1];
campioni_Vd=(trnd(n-3,dim).*sqrt(s2Vd))+Vd; %p=2 ,  gdl = n-p-1 --> n-3
campioni_E=(trnd(n-3,dim).*sqrt(s2E))+emivita_corto; %p=2
campioni_El=(trnd(n-3,dim).*sqrt(s2El))+emivita_lungo; %p=2
campioni_F=(trnd(n-2,dim).*sqrt(s2F))+fraction; %p=1

% SECONDA PARAMETRIZZAZIONE: A1, A2, alpha1 e alpha2 con dose 49650 pmol 

dose=49650;
alpha1=log(2)/emivita_corto; %costante rapida
alpha2=log(2)/emivita_lungo; %costante lenta

% A1+A2=dose/Vd 
% Fraction=A1/(A1+A2)

A1=(fraction*dose)/(Vd*10^3);
A2=(dose/(Vd*10^3)) - A1;

% DISTRIBUZIONI: RICAVATE DA QUELLE DELLA PRIMA PARAMETRIZZAZIONE 

campioni_alpha1=log(2)./campioni_E;
campioni_alpha2=log(2)./campioni_El;
campioni_A1=(campioni_F.*dose)./(campioni_Vd.*10^3);
campioni_A2=(dose./(campioni_Vd.*10^3))-campioni_A1;

% INTERVALLI DI CONFIDENZA BILATERALI CON ALPHA 0.05: RICAVATI CON PERCENTILI

i_alpha1=[prctile(campioni_alpha1,2.5) prctile(campioni_alpha1,97.5)];
i_alpha2=[prctile(campioni_alpha2,2.5) prctile(campioni_alpha2,97.5)];
i_A1=[prctile(campioni_A1,2.5) prctile(campioni_A1,97.5)];
i_A2=[prctile(campioni_A2,2.5) prctile(campioni_A2,97.5)];

% TERZA PARAMETRIZZAZIONE: V1,k01,k12,k21
k12=alpha1-((A1/dose)*Vd*10^3*(alpha1-alpha2));
k01=(alpha1*alpha2)/k12;
k21=alpha1+alpha2-k01-k12;

% DISTRIBUZIONI
campioni_k12=campioni_alpha1-((campioni_A1./dose).*(campioni_Vd.*10^3).*(campioni_alpha1-campioni_alpha2));
campioni_k01=(campioni_alpha1.*campioni_alpha2)./campioni_k12;
campioni_k21=campioni_alpha1+campioni_alpha2-campioni_k01-campioni_k12;

% INTERVALLI DI CONFIDENZA: RICAVATI CON PERCENTILI

i_k12=[prctile(campioni_k12,2.5) prctile(campioni_k12,97.5)];
i_k01=[prctile(campioni_k01,2.5) prctile(campioni_k01,97.5)];
i_k21=[prctile(campioni_k21,2.5) prctile(campioni_k21,97.5)];





  