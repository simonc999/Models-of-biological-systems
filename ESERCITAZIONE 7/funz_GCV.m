function f2 = funz_GCV(gamma,G,D,sigmav,z,n)
uS=(inv(G'*inv(sigmav)*G+gamma*D'*D))*G'*inv(sigmav)*z;
yS=G*uS;
residui2=z-yS;
H=G*inv(G'*inv(sigmav)*G+gamma*D'*D)*G'*inv(sigmav);
q=trace(H);
f2=(n*residui2'*inv(sigmav)*residui2)/((n-q)^2);
end

