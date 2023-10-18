V1 = 3.29;
k01 = 6.54*10^-2; 
k12 = 5.68*10^-2; 
k21 = 7.16*10^-2;
d = 49650*10^(-12); 
d2=(k21-k12)


%matrici della dinamica (compartimento 1)
A1 = [-(k01+k21) k12; k21 -k12];
B1 = [d; 0];
C1 = [1/V1 0];
D=0;
Q1=C1*V1;

%sistema lineare descritto in variabili di stato
sys1 = ss(A1,B1,C1,D);
% sys2 = ss(A2,B2,C2,D2); 

figure(1)
[Y1 T1] = impulse(sys1);
subplot(2,1,1)
plot(T1, Y1)
title('Andamento nel compartimento 1')
xlabel('tempo (min)')
ylabel('concentrazione (mg/l)')
subplot(2,1,2)
semilogy(T1,Y1)
title('Grafico semilogaritmico')
xlabel('tempo (min)')
ylabel('ln(concentrazione)')


%a)Tempo di eliminazione
% tau11=1/(k01-k21+k12)
% tau22=1/(k21-k12)
% Telim11=5*tau11
% Telim22=5*tau22
tau1=1/0.1722
tau2=1/0.0216;
Telim1=5*tau1
Telim2=5*tau2

%b)c(0)=cmax
conc = d/V1;
cmax=max(Y1);


%c)AUC, calcolo area nel compartimento 1
AUC1=trapz(T1,Y1);
AUC2=d/(V1*k01);


%d)Clearance, essendo che il tutto da 1 viene intrappolato nel
%compartimento 2 venendo a creare un loop infinito senza uscita 
Clearance=V1*k01;
AUC3=d/Clearance;
Clearance=d/AUC3;

%Funzione di trasferimento e Analisi del Sistema
trasf=tf(sys1);
p=pole(sys1);
z=zero(sys1);
[R,P,K]=residue(trasf.Numerator{1},trasf.Denominator{1});


