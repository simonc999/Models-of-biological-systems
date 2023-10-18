%% ES 4


clear 
close all


V1 = 5;      %[l]
k01 = 1.2;   %[h^-1]  processi escretori del primo compartimento
k12 = 0.859; %[h^-1]  processi escretori del secondo compartimento
k21 = 2.22;  %[h^-1]
k31 = 0.031; %[h^-1]
k13 = 0.008; %[h^-1] processo di assorbimento
D_es = 500;  %[mg]


%matrici della dinamica (compartimento 1)

A1 = [-(k21+k31+k01), k12, k13; k21 -k12 0; k31 0 -k13];
B1 = [D_es; 0; 0];
C1 = [1/V1 0 0];
D  = [0];

%sistema lineare descritto in variabili di stato
sys1 = ss(A1,B1,C1,D);

[Concentrazione,Tempo] = impulse(sys1);

figure(1)

subplot(3,1,1)
plot(Tempo, Concentrazione)
title('Concentrazione primo compartimento')
xlabel('tempo (h)')
ylabel('concentrazione (mg/l)')
subplot(3,1,2)
semilogy(Tempo, Concentrazione)
title('Concentrazione primo compartimento in scala semilogaritmica')
xlabel('tempo (h)')
ylabel('ln(concentrazione)')
subplot(3,1,3)
plot(Tempo, Concentrazione*V1)
title('Andamento nel primo compartimento')
xlabel('tempo (h)')
ylabel('concentrazione (mg/l)')
%%

AUC1=trapz(Tempo,Concentrazione);

AUMC=trapz(Tempo,Concentrazione.*Tempo);

MRT=AUMC/AUC1;

CL_tot1=D_es/AUC1;

transfer_function = tf(sys1);

figure(100)
subplot(2,1,1)
bode(sys1)
subplot(2,1,2)
pzmap(sys1)

autovalori=eig(A1);

% ABBIAMO UNA QUASI ELIMINAZIONE IN UN POLO, PER QUESTO SI OSSERVANO SOLO
% DUE CONTRIBUTI E NON 3. DOPO AVER ESCLUSO L'AUTOVALORE PIU' BASSO PRENDO
% IL SECONDO PIU BASSO PER CONSIDERARE LA COSTANTE PIU' LENTA 

tau=1/(-autovalori(2));

TEMPO_ELIM=5*tau;

Volume_distribuzione=CL_tot1*tau;

