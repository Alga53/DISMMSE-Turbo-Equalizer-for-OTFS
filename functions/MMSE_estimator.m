function Lext12 = MMSE_estimator(y_DD, H_DD, N0, Tau, Nu, Mean, Var_diag, iter, configs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function realizes our proposed DI-S-MMSE (Turbo) equalizer.
% Input: =========================================== 
%           y_DD ------ the received symbol sequence in the DD domain
%           H_DD  ----- the channel matrix in the DD domain
%           N0 --------- noise spectral density
%           Tau -------- the indices of delay taps
%           Nu --------- the indices of Doppler taps
%           Mean ------ the vector of a priori means of symbols
%           Var_diag --- the diagonal of a priori covariance matrix of symbols 
%           iter --------- the Turbo iteration index
%           configs ----- some predefined parameters]
%
% Output: ==========================================
%           Lext12 ------ the extrinsic LLR passed from the MMSE
%                             estimator to the decoder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lenCode = configs.lenCode;
M = configs.M;
N = configs.N;
tol_gmres = configs.tol_gmres;
Restart = configs.Restart;

lenSym = M*N;
Lext12 = zeros(1, lenCode);          % Initialize the extrinsic information passed from the estimator to the decoder
x_hat = zeros(lenSym, 1);             % Initialize the estimated symbol sequence
xi_seq = zeros(1, lenSym);            % Initialize the vector of xi

% Derive the covariance matrix A = H * V * H' + N_0 * I (refer to equ. 16a)
% Here we use the sparse form of matrix A to save storage and computation
H_DD = sparse(H_DD);
Var = sparse(diag(Var_diag));
idenMtx = sparse( eye(lenSym) );
A = H_DD*Var*H_DD' + N0*idenMtx;

% Apply GMRES at the first outer iteration
if iter == 1
    % Solve the first sparse linear system: A * f_1 = y
    [x1, ~] = gmres(A, y_DD, Restart, tol_gmres, 100);
    
    % Solve the second sparse linear system: A * f_2 = h_n
    h_n = H_DD(:, 1);
    [x2, ~] = gmres(A, h_n, Restart, tol_gmres, 100);
    xi = real(H_DD(:, 1)' * x2);  % xi_n = h_n' * inv(A) * h_n as shown after equ. 19
    xi_seq = xi*ones(1, lenSym);
    % At the first outer iteration, we assume all the xi's have the same value
    % which is validated by Fig. 4.
else
    % Apply FSPAI at the subsequent outer iterations
    A_inv = FSPAI(A, Tau, Nu, configs);
    x1 = A_inv * (y_DD - H_DD * Mean);
    
    for ii = 0:lenSym-1
        % xi_n = h_n' * inv(A) * h_n as shown after equ. 19
        xi_seq(ii+1) = real( H_DD(:, ii+1)' * A_inv * H_DD(:, ii+1) );
    end
end

%%%% Calculate Extrinsic information for Decoder
for nn = 0:lenSym-1
    % Step 1: Derive the estimated symbols (refer to equ. 20)
    xi = xi_seq(nn+1);
    v = Var(nn+1, nn+1);
    x_hat(nn+1) = 1/( 1+(1-v)*xi )*...
        ( H_DD(:, nn+1)' * x1 + xi * Mean(nn+1) );
    
    % Step 2: Calculate the Extrinsic information for the decoder (refer to equ. 21&22)
    Lext12(2*nn+1) = sqrt(8) * ( 1+(1-v)*xi ) * real(x_hat(nn+1)) / (1-v*xi);
    Lext12(2*nn+2) = sqrt(8) * ( 1+(1-v)*xi ) * imag(x_hat(nn+1)) / (1-v*xi);
end

end