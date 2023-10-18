%%%SIMULAZIONE DI MODELLI COMPARTIMENTALI LINEARI - Esercitazione 1

%%%MODELLO A UN COMPARTIMENTO CON ASSORBIMENTO
%parametri
V2 = 5; %litri, volume secondo comparimento
k01 = 1.2; %ore alla meno uno,  processi escretori del primo compartimento
k02 = 1.2; %ore alla meno uno,  processi escretori del secondo compartimento
k21 = 2.2; %ore alla meno uno, processo di assorbimento
d = 500; %mg, bolo per via orale

ka = k21 + k01; %costante di assorbimento

%matrici della dinamica (per vedere andamento nel compartimento 2)
A = [-ka 0; k21 -k02];
B = [d; 0];
C = [0 1/V2];
D = 0;

%matrici della dinamica (per vedere andamento nel compartimento 1)
A1 = [-ka];
B1 = [d];
C1 = [1/V2];

%sistema lineare descritto in variabili di stato
sys = ss(A,B,C,D);
sys1 = ss(A1,B1,C1,D);

figure(1)
[Y1 T1] = impulse(sys1);
subplot(2,1,1)
plot(T1, Y1)
title('Andamento nel primo compartimento')
xlabel('tempo (h)')
ylabel('concentrazione (mg/l)')
subplot(2,1,2)
semilogy(T1,Y1)
title('Grafico semilogaritmico')
xlabel('tempo (h)')
ylabel('ln(concentrazione)');

figure(2)
[Y T] = impulse(sys); %Y= concentrazione T=tempo
subplot(2,1,1)
plot(T,Y)
title('Andamento nel secondo compartimento')
xlabel('tempo (h)')
ylabel('concentrazione (mg/l)')
subplot(2,1,2)
semilogy(T,Y)
title('Grafico semilogaritmico')
xlabel('tempo (h)')
ylabel('ln(concentrazione)')

%(a)Tempo di eliminazione 
tau1=1/(k01+k21);
tau2=1/(k02)
Telim1=5*tau1
Telim2=5*tau2

%Coefficienti valore numeratore e denominatore 
trasf=tf(sys);
[R,P,K]=residue(trasf.Numerator{1},trasf.Denominator{1}); 

%(b)Valore di concentrazione iniziale e concentrazione massima
c0=0;
cmax=max(Y);

%(c) area sotto la curva di concentrazione del compartimento 2
AUC_1 =trapz(T1,Y1);

AUC_2 = trapz(T,Y)
Clearance1=d/AUC_2;
% Clearance2=V2*(k02);
AUC_21=d/Clearance1;
% AUC_22=d/(V2*k02);


%(d)Frazione di farmaco che passa dal compartimento 1 al compartimento 2
%(Circa pari a 64.7%)
F=(k21)/(k21+k01);

%(f)Biodisponibilità
Bio2=AUC_1*(d/2)/d*AUC_2;


%(g)
Velapp=ka;

