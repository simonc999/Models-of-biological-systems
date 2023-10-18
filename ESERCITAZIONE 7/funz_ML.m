function f3 = funz_ML(gamma,G,D,sigmav,z,n)
uS=(inv(G'*inv(sigmav)*G+gamma*D'*D))*G'*inv(sigmav)*z;
yS=G*uS;
residui=z-yS;
H=G*inv(G'*inv(sigmav)*G+gamma*D'*D)*G'*inv(sigmav);
q=trace(H);

c=(residui'*inv(sigmav)*residui*q)/((n-q)*uS'*D'*D*uS);

%Sia a destra che a sinistra ho gamma:FUNZIONE IMPLICITA
f3=abs(gamma-c);
end

