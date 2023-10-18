function dq = ASSORBIMENTO_LIN(q,k01,k02,k21,~,~)

dq=zeros(size(q));

dq(1)=-(k01+k21)*q(1);
dq(2)=k21*q(1)-k02*q(2);
end

