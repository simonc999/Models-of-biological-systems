function dp = INCOMPETITIVA(p, k1in, k_1in, k2in, k3in, k_3in)

dp = zeros(size(p));

dp(1)= -k1in*p(2)*p(1) +k_1in*p(3);                                            %1 s
dp(2)= -k1in*p(2)*p(1) +(k_1in+k2in)*p(3);                                     %2 e
dp(3)= k1in*p(2)*p(1) -(k_1in+k2in)*p(3) +k_3in*p(6) -k3in*p(3)*p(5);          %3 c  (ES)
dp(4)= k2in*p(3);                                                              %4 p
dp(5)= -k3in*p(3)*p(5) +k_3in*p(6);                                            %5 i
dp(6)= k3in*p(3)*p(5) -k_3in*p(6);                                             %6 p2 (EIS)

end
