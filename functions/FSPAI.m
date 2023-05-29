function A_inv = FSPAI(A, Tau, Nu, configs)
M = configs.M;
N = configs.N;
tol_A = configs.tol_A;
maxDeg = configs.maxDeg;
tol_fspai = configs.tol_fspai;
maxiter_fspai = configs.maxiter_fspai;

P = length(Tau);
diag_A_inv = sparse( diag( 1./sqrt( diag(A) ) ) );
A_norm = diag_A_inv*A*diag_A_inv;

Edge_A = zeros(M*N, P*(P-1)+1);
for kk = 0:M*N-1
    m = mod(kk, M);     % delay
    n = floor(kk/M);      % Doppler
    
    pt1 = 1; 
    for pp1 = 0:P-1
        tau1 = Tau(pp1+1);
        nu1 = Nu(pp1+1);
        
        for pp2 = 0:P-1
            if pp1 == 0 || (pp1 > 0 && pp2 ~= pp1)
                tau2 = Tau(pp2+1);
                nu2= Nu(pp2+1);
                
                delta_tau = tau1-tau2;
                delta_nu = nu1-nu2;
                
                posi_x = mod(m+delta_tau, M)+M*mod(n+delta_nu, N)+1;
                
                if abs(A_norm(posi_x, kk+1)) > tol_A
                    if min( abs (Edge_A(kk+1, :)-posi_x) ) ~= 0
                        Edge_A(kk+1, pt1) =  posi_x;
                        pt1 = pt1+1;
                    end
                end
            end
        end
    end
    
    if pt1 <= maxDeg+1
        Edge_A(kk+1, 2:pt1) = zeros(1, pt1-1);
    end
end

%%%%%%%% FSPAI
L = eye(M*N);

for nn = 0:M*N-1
    idx_J = nn+1;
    idx_J_tilde = []; % The set of indices except 'k'
    
    a = real( A(nn+1, nn+1) );
    
    L(nn+1, nn+1) = 1 / sqrt( a );
    len = 0;
    while (len < maxiter_fspai)
        
        idx_J_hat = zeros(1, P^4);
        len_IdxJhat = 0;
        for ee = 1:length(idx_J)
            pt = 1;
            
            idx = Edge_A(idx_J(ee), pt);
            while ( idx ~= 0 && pt < P*(P-1)+1  )
                if idx > nn+1 && isempty( find(idx_J==idx, 1) ) &&...
                        isempty( find( idx_J_hat(1:len_IdxJhat)==idx, 1 ) )
                    idx_J_hat(len_IdxJhat+1) = idx;
                    len_IdxJhat = len_IdxJhat+1;
                end
                pt = pt+1;
                idx = Edge_A(idx_J(ee), pt);
            end
            
        end
        idx_J_hat = idx_J_hat(1:len_IdxJhat);
        
        eta = zeros(1, len_IdxJhat);
        for jj = 1:len_IdxJhat
            eta(jj) = abs( A(idx_J, idx_J_hat(jj))' * L(idx_J, nn+1) )^2 / a;
        end
        
        if isempty(eta)
            break;
        else
            [val_max, idx_max] = max( eta );
            if max(val_max) <= tol_fspai
                break;
            end
        end
        
        idx_J = [idx_J, idx_J_hat(idx_max)];
        idx_J_tilde = idx_J(2:end);
        
        yk = A(idx_J_tilde, idx_J_tilde) \ A(idx_J_tilde, nn+1);
        L(nn+1, nn+1) = 1 / sqrt( a-A(idx_J_tilde, nn+1)'*yk );
        L(idx_J_tilde, nn+1) = -L(nn+1, nn+1)*yk;
        len = length(idx_J_tilde);
    end
end

L1 = sparse(L);
L2 = sparse(L');
A_inv = L1*L2;

end