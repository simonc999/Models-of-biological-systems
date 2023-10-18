function f1 = funz_discrepanza(gamma,G,D,sigmav,z,traccia)

uS=(inv(G'*inv(sigmav)*G+gamma*D'*D))*G'*inv(sigmav)*z;
yS=G*uS;
residui2=z-yS;
RSS=sum(residui2.^2);

f1=abs(RSS-traccia);

end

