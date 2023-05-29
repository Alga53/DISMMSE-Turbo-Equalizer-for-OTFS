function x_T = ISFFT_Heisenberg(x_DD, M, N)

R = reshape(x_DD, M, N);
for mm = 0:M-1
    R(mm+1, :) = sqrt(N) * ifft( R(mm+1, :) );
end
x_T = reshape(R, M*N, 1);

end