function dq = ASSORBIMENTO_NO_LIN_MM(q,k01,k02,km,Vmax)

dq=zeros(size(q));

k21=Vmax/(km+q(1));

dq(1)=-(k01+k21)*q(1);
dq(2)=k21*q(1)-k02*q(2);

end

