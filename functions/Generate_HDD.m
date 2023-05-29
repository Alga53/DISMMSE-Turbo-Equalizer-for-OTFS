function H_DD = Generate_HDD(M, N, h, Tau, Nu)

P = length(h);
Edge_y = zeros(M*N, P);
Edge_x = zeros(M*N, P);

for pp = 0:P-1
    tau = Tau(pp+1);
    nu = Nu(pp+1);
    
    for kk = 0:M*N-1
        m = mod(kk, M);     % delay
        n = floor(kk/M);      % Doppler
        Edge_y(kk+1, pp+1) =  mod(m-tau, M)+M*mod(n-nu, N)+1;
        Edge_x(kk+1, pp+1) =  mod(m+tau, M)+M*mod(n+nu, N)+1;
    end
end

H_DD = zeros(M*N, M*N);

h_rotat = zeros(1, P);
for pp = 0:P-1
    h_rotat(pp+1) = h(pp+1)*exp(-1i*2*pi*Tau(pp+1)*Nu(pp+1)/M/N);
end
pt = 1;
for rr = 0:M*N-1
    n = floor(rr/M);
    m = rr-n*M;
    for pp = 0:P-1
        tau = Tau(pp+1);
        nu = Nu(pp+1);
        
        m_tau = mod(m-tau, M);
        
        H_DD(rr+1, Edge_y(rr+1, pp+1)) = h(pp+1)*...
            exp(1i*2*pi/M/N*nu*(m_tau));
        
        if m < tau
            H_DD(rr+1, Edge_y(rr+1, pp+1)) = H_DD(rr+1, Edge_y(rr+1, pp+1))*...
                exp(-1i*2*pi*n/N);
        end
    end
end

end