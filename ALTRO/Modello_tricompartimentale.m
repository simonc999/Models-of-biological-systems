V2 = 5; %litri, volume secondo comparimento
k01 = 1.2; %ore alla meno uno,  processi escretori del primo compartimento
k12 = 0.859; %ore alla meno uno,  processi escretori del secondo compartimento
k21 = 2.22;
k31 = 0.031;
k13 = 0.008;%ore alla meno uno, processo di assorbimento
d = 500; %mg, bolo per via orale


%matrici della dinamica (compartimento 1)
A1 = [-(k21+k31+k01), k12, k13; k21 -k12 0; k31 0 -k13];
B1 = [d; 0; 0];
C1 = [1/V2 0 0];
D=[0];


%sistema lineare descritto in variabili di stato
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
ylabel('ln(concentrazione)')


AUC1=trapz(T1,Y1);
AUMC=trapz(T1,Y1.*T1);
MRT=AUMC/AUC1
CL_tot1=d/AUC1
vd1=CL_tot1/4.0489
trasf=tf(sys1);
p=pole(sys1);
z=zero(sys1);
tau1=1/4.0489