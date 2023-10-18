% la funzione definisce la CINETICA ENZIMATICA a DUE STADI [p.29]
% x = [s e c p]

function dc = ENZIMA(x,k1e,k2e,k_1e,k_2e)

dc = zeros(size(x));


dc(1) = -k1e*x(2)*x(1)+k_1e*x(3);                        % ds
dc(2) = -k1e*x(2)*x(1)+(k_1e+k2e)*x(3)-k_2e*x(2)*x(4);   % de
dc(3) = k1e*x(2)*x(1)-(k_1e+k2e)*x(3)+k_2e*x(2)*x(4);    % dc
dc(4) = k2e*x(3)-k_2e*x(2)*x(4);                         % dp

end
