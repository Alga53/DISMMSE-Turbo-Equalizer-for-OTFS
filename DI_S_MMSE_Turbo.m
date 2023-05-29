
clc;
close all;
addpath('functions')

%%%%%%%% = Parameters Initialization = %%%%%%%%%%
%%% Parameters of transmitter
M = 32;                % the number of subcarriers
N = 16;                 % the number of time slots
lenCP = 16;            % the length of CP per OTFS frame, lenCP > tau_max
P = 4;                   % the number of reflectors

% Encoder
Rc = 1/2;               % Code Rate
mem = 2;              % Memory order of the convolutional encoder

% Modulator
bitSet = [0 0; 1 0; 0 1; 1 1];
symSet = 1/sqrt(2)*[1+1i, -1+1i, 1-1i, -1-1i];  % Symbol Alphabet
numConste = 4;                                    % the number of constellations
order = log2(numConste);
lenBit = M*N*order*Rc-mem;                  % the length of bit sequence
lenBit_tailed = M*N*order*Rc;                 % the length of bit sequence after adding trellis termination bits
lenCode = M*N*order;                           % the length of code sequence
lenSym = M*N;                                     % the length of symbol sequence

% Interleaver
row_Intlvr = 64;   % Parameters of the rectangular interleaver
col_intlvr = lenCode/row_Intlvr;

%%% Parameters of channel
tau_max = 10;              % the maximum delay, tau_max <= M-1
nu_max = 6;                % the maximum Doppler, nu_max <= N-1

%%% Parameters of receiver
% Sparsification
tol_A = 1e-3;               % the threshold of Sparsification Guideline 1
maxDeg = P/4;            % the threshold of Sparsification Guideline 2
% GMRES
tol_gmres = 1e-3;         % the drop tolerance of GMRES
Restart = 5;                 % the restart parameter of GMRES
% FSPAI
tol_fspai = 1e-3;            % the drop tolerance of FSPAI
maxiter_fspai = P;          % the maximum node of degree of FSPAI

%%% Parameters of simulation
EbN0_dB = 8:9;
EbN0 = 10.^(EbN0_dB/10);
Es = 1;                                                  % the average energy of symbols
Eb = (lenSym+lenCP)*Es/lenBit;                % the average energy of bits
sumSim = 500;                                        % the number of simulations at each SNR
iterTimes = 5;                                        % the number of Turbo iterations
BER = zeros(iterTimes, length(EbN0));

% Create a stucture array to store neccessary parameters 
% to facilitate passing arguments 
configs.M = M;
configs.N = N;
configs.lenCode = lenCode;
configs.tol_A = tol_A;
configs.maxDeg = maxDeg;
configs.tol_gmres = tol_gmres;
configs.Restart = Restart;
configs.tol_fspai = tol_fspai;
configs.maxiter_fspai = maxiter_fspai;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for snr = 1:length(EbN0_dB)
    N0 = Eb/EbN0(snr);
    
    %%% If you would like to speed up this code, you may use "parfor sim = 1:sumSim" 
    %%% so as to execute for-loop on a parallel pool of workers on your
    %%% multi-core computer.
    for sim = 1:sumSim
        %%%%%%%% = Transmitter = %%%%%%%%%%
        % = Generating Bit Stream = %
        bit_seq = 1/2+1/2*sign( randn(lenBit_tailed, 1) );
        bit_seq(lenBit+1:end) = zeros(1, mem);      
        % Trellis termination: force the final state of encoder to the all-zero state,
        % which is leveraged in the backward recursion of BCJR algorithm.
        % This overhead caused by trellis termination is also considered. 
        
        % = CC Encoder = %
        code_seq = CC57(bit_seq);
        
        % = Interleaver = %
        code_intlvr = Rect_Interlvr( code_seq, row_Intlvr, col_intlvr);
        
        % = Modulation = %
        x_DD = zeros(lenSym, 1);     % The DD domain transmitted symbol sequence
        for nn = 0:lenSym-1
            [~, posi] = min( sum( abs(code_intlvr(2*nn+1:2*nn+2)'-bitSet), 2) );
            x_DD(nn+1) = symSet(posi);
        end  
        
        % = ISFFT & Heisenberg Transform = %
        x_T = ISFFT_Heisenberg(x_DD, M, N);   % The time domain transmitted symbol sequence
        
        
        %%%%%%%%% = Channel = %%%%%%%%%%%
        % Randomly generate the channel gain, delay and Doppler shift
        [h, Tau, Nu] = CSI_Generator(P, tau_max, nu_max);
        % Derive the channel matrix in the time domain (ref equ. 7)
        H_T = Generate_HT(M, N, h, Tau, Nu);
        % Derive the channel matrix in the DD domain (ref equ. 13)
        H_DD = Generate_HDD(M, N, h, Tau, Nu);
        % Generate AWGN
        n_T = sqrt(N0/2) * (randn(M*N, 1)+1i*randn(M*N, 1));
        
        
        %%%%%%%%% = Receiver = %%%%%%%%%%%
        y_T = H_T*x_T+n_T;
        
        % = Wigner Transform & SFFT = %
        y_DD = Wigner_SFFT(y_T, M, N);
        
        % = DI-S-MMSE Equalizer = %
        Bit_decod = zeros(lenBit, 1);      % Initialize the decoded bit sequence
        Mean = zeros(M*N, 1);             % Initialize the vector of means of the transmitted symbols
        Var = diag( ones(1, M*N) );       % Initialize the covariance matrix of the transmitted symbols
        numErrors = zeros(iterTimes, length(EbN0)); % Initialize the BER of the current simulation
        for iter = 1:iterTimes
            % = MMSE Estimator = %
            Lext12 = MMSE_estimator(y_DD, H_DD, N0, Tau, Nu, Mean, diag(Var), iter, configs);
            
            % = Deinterleaver = %
            Lext12_deintlvr = Rect_Deinterlvr(Lext12, row_Intlvr, col_intlvr);
            
            % = Decoder using BCJR Algorithm = %
            [Lext21_deintlvr, Lapp] = MAPDecoder(Lext12_deintlvr);
            
            % = Interleaver = %
            Lext21 = Rect_Interlvr( Lext21_deintlvr, row_Intlvr, col_intlvr);
            
            % Step 3: Update the means and variances of each symbol
            for nn = 0:M*N-1
                Mean(nn+1) = 1/sqrt(2)*( tanh(Lext21(2*nn+1)/2)+1i*tanh(Lext21(2*nn+2)/2) );
                Var(nn+1, nn+1) = 1-abs(Mean(nn+1))^2;
            end
            
            % = BER calculation = %
            for nn = 0:lenBit-1
                Bit_decod(nn+1) = 1/2*sign(-Lapp(nn+1))+1/2;
            end
            error = sum(abs(Bit_decod-bit_seq(1:lenBit)));
            
            %%%% BER calculation
            numErrors(iter, snr) = error;
            
        end
        BER = BER + numErrors/lenBit/sumSim;
        
        %%% If you use parfor to speed up the calculation, please comment out
        %%% the following code, which may cause errors in parfor loop. 
        clc
        disp('===========================================================')
        display(EbN0_dB, 'EbN0 (dB)');
        display(BER(1, :),'BER the 1st iteration');
        display(BER(2, :),'BER of the 2nd iteration');
        display(BER(5, :),'BER of the 5th iteration');
        disp('===========================================================')
        
        %%% If you use parfor to speed up the calculation, you may
        %%% uncomment the following code to display current EbN0 and
        %%% simulation index. This will allow you to monitor the progress of the simulation.
%         clc
%         display(EbN0_dB(snr), 'Current EbN0 (dB)');
%         display(sim, 'Current simulation index');
%         disp('===========================================================')
        
    end
    
end

figure(1)
semilogy(EbN0_dB,BER(1,:),'-o','LineWidth',2,'Color',[0.25 0.41 0.88]);
hold on;
semilogy(EbN0_dB,BER(2,:),'-o','LineWidth',2,'Color',[0.24 0.57 0.25]);
hold on;
semilogy(EbN0_dB,BER(3,:),'-o','LineWidth',2,'Color',[1 0.5 0.31]);
hold on;
semilogy(EbN0_dB,BER(iterTimes,:),'-o','LineWidth',2,'Color',[0 0 0]);
hold on;
grid on;
str = ['iter','=', num2str(iterTimes)];
legend('iter=1', 'iter=2', 'iter=3', str);
axis([0,15,10^(-5),1]);
set(gcf, 'Color', [1,1,1]);
set(gca, 'Fontname', 'Times New Roman','FontSize',13);
xlabel('$E_{b} / N_{0}\ (\mathrm{dB})$','interpreter','latex','fontsize',14);
ylabel('BER','fontsize',14);