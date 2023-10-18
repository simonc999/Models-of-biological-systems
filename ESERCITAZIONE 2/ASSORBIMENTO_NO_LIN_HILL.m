function dq = ASSORBIMENTO_NO_LIN_HILL(q,k01,k02,km,Vmax,Q)

dq=zeros(size(q));

k21=(Vmax*(q(1)^(Q-1)))/((km^(Q))+(q(1)^(Q)));

dq(1)=-(k01+k21)*q(1);
dq(2)=k21*q(1)-k02*q(2);

end

