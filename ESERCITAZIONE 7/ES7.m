clear 
close all
clc
% OBBIETTIVO: VALUTARE SECREZIONE DELL'INSULINA IN RISPOSTA AD UNO STIMOLO
% DI GLUCOSIO (VALUTAZIONE EFFETTUATA PER RICERCARE PAZIENTI DIABETICI CHE
% NON REAGISCE MOLTO CON LA PRODUZIONE DI INSULINA A SEGUITO DI STIMOLI)

% CON LA STIMOLAZIONE CON IL GLUCOSIO CAMBIO LA PRODUZIONE DELLA NOSTRA
% INSULINA E MISURANDONE L'USCITA CON LA DECONVOLUZIONE VOGLIO RICOSTRUIRNE
% L'INGRESSO

% INGRESSO --> INFUSIONE ENDOVENOSA UP-DOWN
% AUMENTO INIZIALMENTE E POI DIMINUISCO
% POI VALUTO CONC PLASMATICA DI C-PEPTIDE (Y) E TRAMITE DECONVOLUZIONE
% RICOSTRUISCO U 

% PUO ESSERE APPLICATO SOLO AI SISTEMI DINAMICI LINEARI, PERCHE' HO
% USCITA CHE E' INTEGRALE DI CONVOLUZIONE DELLA RISPOSTA IMPULSIVA CON
% L'INGRESSO SOLO SE HO SIS DINAMICO LINEARE E NON HO LA MATRICE DELLA
% DINAMICA D

yb= 0.32;

% CARICO I DATI

datiEs7;              % Dati
stimaEs7;             % Parametri soggetto 
t = DATA(:,1);        % Tempi [min]
c = DATA(:,2);        % Concentrazioni [pmol/ml]

%% %%%%%%%%%%%%%%%%%%%%%%% RAW DECONVOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  
% y = G*u

% COSTRUISCO MATRICE G (INTEGRALE DI CONVOLUZIONE)

% N = NUMERO DI MISURE

% M = DISCRETIZZAZIONE (GRIGLIA VIRTUALE)

alfa = log(2)/temp_emi_short;
beta = log(2)/temp_emi_long;
A = fraction/vol;
B = (1/vol)-A; 

n = length(c);
t_comp=[0;t];

g = zeros(1,length(t));
figure(1)


%         COSTRUZIONE MATRICE G E INTEGRALE DI CONVOLUZIONE

for i = 1:length(t_comp)-1
    int = [t_comp(i):10:t_comp(i+1)];   % ESTREMI E PASSO DI INTEGRAZIONE

    g=A*exp(-int*alfa)+B*exp(-int*beta);

    G_1(i) = trapz(int,g); % CALCOLO L'INTEGRALE DI CONVOLUZIONE
   
    plot(int,g)
    xlabel('tempo [min]');
    ylabel('g(t) [1/ml]');
    title('Risposta Impulsiva');
    hold on
end

% CALCOLO LA MATRICE G MATRICE DI CONVOLUZIONE, TRIANGOLARE INFERIORE

for i = 1:n
    for j = i:n
        G(j,i) = G_1(j-i+1);
    end
end

% METODO ALLE DIFFERENZE

% TOLGO IL BASALE

z = c - yb;

% RICOSTRUZIONE INGRESSO

uR=G^-1*z;

% RICONVOLUZIONE

z_riconv = G*uR;

% CALCOLO I RESIDUI

res_raw = z - z_riconv;

%CALCOLO DEI RESIDUI CON IL MODELLO A CV COSTANTE 

CV=0.04;

x2=  [20; 1.2; 10; 0.5];
%options = optimoptions('lsqnonlin','Display','iter');

sigma_k=CV.*c;

res_pes=res_raw./sigma_k;

figure(2)
subplot(4,1,1)
stairs(t,uR)
hold on
stairs(t,uR,'o')
xlabel("Tempo (min)")
ylabel("Ingresso (pmol)")
title("Raw Convolution")
grid on; 

subplot(4,1,2)
plot(t,z_riconv,t,z,'o') 
xlabel("Tempo (min)")
ylabel("uscita (pmol/ml)")
title("Riconvoluzione e dati originari")
grid on; 

subplot(4,1,3)
plot(t,res_raw,t,zeros(length(t),1),'o') 
xlabel("Tempo (min)")
ylabel("Residui(pmol/ml)")
title("Residui"), 
ylim([-0.001,0.001])
grid on;

subplot(4,1,4)
plot(t,res_pes,t,zeros(length(t),1),'o')
xlabel("Tempo (min)")
ylabel("Residui(pmol/ml)")
title("Residui Modello a CV Costante"), 
ylim([-0.001,0.001]), grid on;

%% %%%%%%%%%%%%%%%%%%%%%%% REGOLARIZZAZIONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VOGLIO TENER CONTO DEL RUMORE NEI DATI, E PENALIZZARE LE SOLUZIONI POCO
% REGOLARI. DEVO INTRODURRE MODELLO DELL'ERRORE, E MATRICE DI PENALITÀ.
% POSSO INOLTRE INFITTIRE LA GRIGLIA VIRTUALE SU CUI RICOSTRUIRE IL SEGNALE


griglia_virtuale = [ 0 : 1 : max(t) ]; %max(t) è 240 
griglia_virtuale =griglia_virtuale';


figure(3)
for i = 1:length(griglia_virtuale)-1
    %Il passo viene calibrato a seconda della preicisione desiderata
    int = [griglia_virtuale(i):0.01:griglia_virtuale(i+1)];
    g_reg(i) = trapz(int,A*exp(-int*alfa)+B*exp(-int*beta));
    plot(int, A*exp(-int*alfa)+B*exp(-int*beta))
    xlabel('tempo [min]');
    ylabel('g(t) [1/ml]');
    title('Risposta Impulsiva');
    hold on
end

%PASSO 1:Costruisco come prima G, inizialmente quadrata

for i = 1:length(griglia_virtuale)-1
    for j = i:length(griglia_virtuale)-1
        G_reg(j,i) = g_reg(j-i+1);
    end
end
 
% PASSO 2:Considero esclusivamente gli istanti di campionamento che  mi
% corrispondono all'uscita (0,10,20,30), eliminando tutti gli altri.

GR = G_reg(t_comp(2:end),:);

% Costruisco la MATRICE DI PENALITA', limitandomi allo studio della 
% derivata prima.

% P*ù SIA IL VETTORE DELLE DIFFERENZE KAPPESIME DI ù QUINDI CIRCA KAPPESIMA
% DERIVATA DISCRETA DI ù

d = [ 1;-1;zeros(length(griglia_virtuale)-3,1)];
D = tril(toeplitz(d));
sigmav = diag((CV.*c).^2);   % MANTENGO SEMPRE IL BASALE PER ERRORE DI MISURA
traccia = trace(sigmav);

% DERIVATA ALLE DIFFERENZE 

K_raw=sqrt(max(eig(G'*G)/min(eig(G'*G))));

% PIU IMPONGO LE DERIVATE PICCOLE SALENDO DI GRADO DERIVATIVO PIU IMPONGO
% LA SALITA LENTAMENTE

%% %%%%%%%%%%%%%%%%%%%%%%% SCELTA DEL GAMMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%% %%%%%%%%%%%%%%%%%%%%% criterio di discrepanza %%%%%%%%%%%%%%%%%%%%%%%%%%

% IL CRITERIO DI DISCREPANZA DI TWOMEY È INTUITIVO MA SI PUÒ 
% DIMOSTRARE CHE, IN MEDIA, CONDUCE AD OVERSMOOTHING (Γ TROPPO 
% GRANDE).

% OBBIETTIVO: SODDISFARE UGUAGLIANZA RSS = TRACCIA SIGMAV

%  Scegliere γ

%  calcolare la stima di u per il γ scelto

%  Calcolare la riconvoluzione

%  Calcolare il vettore dei residui

%  Calcolare RSS

%  fino a quando soddisfo circa l’uguaglianza
%  parto con alcuni valori di gamma in scala log
 

vett_gamma=logspace(-10,6,100);   %  Scegliere γ
for i=1:length(vett_gamma)

    % DA FORMULA RICALCOLO LA STIMA DI u
    
    %  calcolare la stima di u per il γ scelto
    ui=(inv(GR'*inv(sigmav)*GR+vett_gamma(i)*D'*D))*GR'*inv(sigmav)*z;

    %  Calcolare la riconvoluzione
    yi=GR*ui;

    %  Calcolare il vettore dei residui
    residui_i=z-yi;

    %  Calcolare RSS
    RSS_D=sum(residui_i.^2);
    
    % SODDISFARE L'UGUAGLIANZA EQUIVALE A MINIMIZZARE LA DIFFERENZA
    fD(i)=RSS_D-traccia;
end

[diff_minore,idx_gammaOtt]=min(abs(fD));

gammaOtt=vett_gamma(idx_gammaOtt);

%partendo dal valore minimo trovato con il ciclo precedente, provo altri
%valori di gamma rimanendo intorno a questo valore

vett_gamma2=linspace(vett_gamma(idx_gammaOtt-2),vett_gamma(idx_gammaOtt+2),100);

for i=1:length(vett_gamma2)
    ui=(inv(GR'*inv(sigmav)*GR+vett_gamma2(i)*D'*D))*GR'*inv(sigmav)*z;
    yi=GR*ui;
    residui_i=z-yi;
    RSS_D=sum(residui_i.^2);
    fD2(i)=RSS_D-traccia;
end

[diff_minore2,idx_gammaOtt2]=min(abs(fD2));
gammaOtt2=vett_gamma2(idx_gammaOtt2);

figure(4)

subplot(1,2,1);            % LOGSPACE
plot(fD);
grid on;
ylabel('|RSS-traccia|');
xlabel('indici');
title('scala log');

subplot(1,2,2);            % LINSPACE
plot(fD2);
grid on;
ylabel('|RSS-traccia|');
xlabel('indici');
title('scala lin');
sgtitle('|RSS - traccia(sigmav)| per diversi valori di gamma');


%come verifica utilizzo un metodo di minimizzazione

gamma0=0.005; 

%parto da un valore vicino al gamma ottimo trovato prima


gammaD=fminunc(@(gammaD) funz_discrepanza(gammaD,GR,D,sigmav,z,traccia),gamma0);
uS=(inv(GR'*inv(sigmav)*GR+gammaD*D'*D))*GR'*inv(sigmav)*z;
yS=GR*uS;
residui2=z-yS;
dev_st=sqrt(sigmav);
residui2_pes=residui2./dev_st;
dev_stResidui=std(residui2_pes);

figure(5)
subplot(3,1,1);                                    % INGRESSO DISCREPANZA
stairs(griglia_virtuale(2:end,:),uS);
grid on;
xlabel('tempo [min]');
ylabel('u [pmol]');
title(['Regolarizzazione, gamma = ' num2str(gammaD)]);

subplot(3,1,2);                                    % USCITA DISCREPANZA
plot(t,yS,t,z,'*');
grid on;
xlabel('tempo [min]');
ylabel('concentrazione [pmol/ml]');
title('Riconvoluzione');

subplot(3,1,3);                                    % RESIDUI PESATI DISCREPANZA
plot(t,residui2_pes,'*','MarkerEdgeColor','b');
grid on;
xlabel('tempo [min]');
ylabel('res');
title('Residui Pesati');
sgtitle('Criterio Discrepanza');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% criterio GCV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% UTILIZZO SOLO fminunc

gammaGCV=fminunc(@(gammaGCV) funz_GCV(gammaGCV,GR,D,sigmav,z,n),gamma0);
uS2=(inv(GR'*inv(sigmav)*GR+gammaGCV*D'*D))*GR'*inv(sigmav)*z;
yS2=GR*uS2;
residui3=z-yS2;
residui3_pes=residui3./dev_st;
 
figure(6);
subplot(3,1,1);                                    % INGRESSO GCV
stairs(griglia_virtuale(2:end,:),uS2);
grid on;
xlabel('tempo [min]');
ylabel('u [pmol]');
title(['Regolarizzazione, gamma = ' num2str(gammaGCV)]);


subplot(3,1,2);                                    % USCITA GCV
plot(t,yS2,t,z,'*');
grid on;
xlabel('tempo [min]');
ylabel('concentrazione [pmol/ml]');
title('Riconvoluzione');


subplot(3,1,3);                                    % RESIDUI GCV
plot(t,residui3_pes,'*','MarkerEdgeColor','b');
hold on
yline(0)
grid on;
xlabel('tempo [min]');
ylabel('res');
title('Residui Pesati');
sgtitle('Criterio GCV');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% criterio ML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UTILIZZO SOLO fminunc

gammaML=fminunc(@(gammaML) funz_ML(gammaML,GR,D,sigmav,z,n),gamma0);
uS3=(inv(GR'*inv(sigmav)*GR+gammaML*D'*D))*GR'*inv(sigmav)*z;
yS3=GR*uS3;
residui4=z-yS3;
residui4_pes=residui4./dev_st;
 
figure(7);

subplot(3,1,1);                                    % INGRESSO ML
stairs(griglia_virtuale(2:end,:),uS3);
grid on;
xlabel('tempo [min]');
ylabel('u [pmol]');
title(['Regolarizzazione, gamma = ' num2str(gammaML)]);


subplot(3,1,2);                                    % USCITA ML
plot(t,yS3,t,z,'*');
grid on;
xlabel('tempo [min]');
ylabel('concentrazione [pmol/ml]');
title('Riconvoluzione');


subplot(3,1,3);                                    % RESIDUI ML
plot(t,residui4_pes','*','MarkerEdgeColor','b');
hold on
yline(0)
grid on;
xlabel('tempo [min]');
ylabel('res');
title('Residui Pesati');
sgtitle('Criterio ML');
