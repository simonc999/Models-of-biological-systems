% Eta': 32
% Stato salute: Normale
% Sesso: F
% Altezza: 1.64
% Peso: 59.6

% Parametri stimati sul soggetto con dose unitaria essendo che si desidera
% ricostruire un impulso unitario g(t)

beta_vd = [1.0912,1.8623,-0.0101]; % VEDI BETA GENERALE VOLUME DI DISTRIBUZIONE

beta_t_emi_s = [5.8970,-0.0280,-0.0088]; % VEDI BETA GENERALE T EMI VELOCE

beta_t_emi_l = [44.2966,0.1696,-9.4406]; % VEDI BETA GENERALE T EMI LENTO

beta_fraction = [0.7401,0.0011]; % VEDI BETA GENERALE FRACTION


vol = 3.8229*10^(3);            %Vedi Vdistr
temp_emi_short = 5.1746; %Vedi Tcortol
temp_emi_long = 32.1870; %Vedi Tlungol
fraction = 0.7617;      %Vedi Fl 
 


