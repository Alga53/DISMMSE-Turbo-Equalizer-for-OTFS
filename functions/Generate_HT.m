function H_T = Generate_HT(M, N, h, Tau, Nu)
P = length(Tau);
H_T = zeros(M*N, M*N);

Pai = eye(M*N);
Delta = zeros(1, M*N);
for ii = 0:M*N-1
    Delta(1, ii+1) = exp(1i*2*pi*ii/M/N);
end

for pp = 0:P-1
    Pai_pp = circshift(Pai, -Tau(pp+1), 2);
    Delta_pp = diag( Delta.^Nu(pp+1) );
    
    H_T = H_T + h(pp+1) * Pai_pp * Delta_pp;
end

end