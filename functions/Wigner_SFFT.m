function y_DD = Wigner_SFFT(y_T, M, N)

R = reshape(y_T, M, N);
for mm = 0:M-1
    R(mm+1, :) = sqrt(1/N) * fft( R(mm+1, :) );
end
y_DD = reshape(R, M*N, 1);

end