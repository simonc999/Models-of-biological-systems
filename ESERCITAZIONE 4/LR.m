% LR sta per Linear Residual

function [RES CONC SYS]= LR (par,t,z,pesi,dose)

k01 = par(1);
k12 = par(2);
k21 = par(3);
V = par(4);

% definisco le matrici e creo il sistema 
A = [(-k01-k21) k12;k21 -k12];
B = [dose;0];
C = [1/V 0];
D = [0];

SYS = ss(A,B,C,D);

[c_imp, t_imp] = impulse(SYS,1:t(end));

% per calcolare i residui ho tempi differenti, perciò uso la funzione
% 'interp1' che mi permette di interpolare i valori correttamente

CONC = interp1(t_imp,c_imp,t);
RES = (pesi)*(z-CONC);

end 