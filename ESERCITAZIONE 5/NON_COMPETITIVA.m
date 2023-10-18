function dp = NON_COMPETITIVA(p, k1in, k_1in, k2in, k_2in, k3in, k_3in, k4in, k_4in)

dp = zeros(size(p));

dp(1)= -k1in*p(2)*p(1) +k_1in*p(3);                            %1 s
dp(2)= -k1in*p(2)*p(1) -k3in*p(2)*p(5) +(k_1in+k2in)*p(3)  ... %2 e
       +(k_3in+k4in)*p(6) -k_2in*p(2)*p(4) -k_4in*p(2)*p(7); 
dp(3)= k1in*p(2)*p(1) -(k_1in+k2in)*p(3) +k_2in*p(2)*p(4);     %3 c  (ES)
dp(4)= k2in*p(3) -k_2in*p(2)*p(4);                             %4 p
dp(5)= -k3in*p(2)*p(5) +k_3in*p(6);                            %5 i
dp(6)= k3in*p(2)*p(5) -(k_3in+k4in)*p(6) +k_4in*p(2)*p(7);     %6 c2 (EI)
dp(7)= k4in*p(6) -k_4in*p(2)*p(7);                             %7 p2 (EIS)

end
