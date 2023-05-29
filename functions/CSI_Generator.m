function [h, Tau, Nu] = CSI_Generator(P, tau_max, nu_max)

h = 1/sqrt( 2*P )*( randn(1, P) + 1i*randn(1, P) );
Reflec_idx = randperm( (2 * nu_max + 1) * (tau_max + 1), P) - 1;
Tau = mod(Reflec_idx, tau_max+1);
Nu = floor(Reflec_idx ./ (tau_max+1)) - nu_max;

end