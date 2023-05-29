function   [Lext21_deintlvr, Lapp] = MAPDecoder(Lext12_deintlvr)

lenBit = length(Lext12_deintlvr)/2;
Lext21_deintlvr = zeros(1, 2*lenBit);
Lapp = zeros(1, lenBit);
P0 = 1/2 * ( 1+tanh(Lext12_deintlvr/2) );
P1 = 1-P0;

cProb = zeros(4, 4, lenBit);  % The initialization of Gamma
aProb = zeros(4, lenBit+1);  % The initialization of Alpha
bProb = zeros(4, lenBit+1);  % The initialization of Beta
cProb_b1 = zeros(4, 4);
cProb_b2 = zeros(4, 4);

CnctMtx = [1 0 1 0; 1 0 1 0; 0 1 0 1; 0 1 0 1];
CnctMtx_a_0 = [1 0 0 0; 1 0 0 0; 0 1 0 0; 0 1 0 0];
CnctMtx_a_1 = CnctMtx-CnctMtx_a_0;
CnctMtx_b1_0 = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
CnctMtx_b1_1 = CnctMtx-CnctMtx_b1_0;
CnctMtx_b2_0 = [1 0 0 0; 0 0 1 0; 0 0 0 1; 0 1 0 0];
CnctMtx_b2_1 = CnctMtx-CnctMtx_b2_0;

% Calculation of Gamma
for kk = 0:lenBit-1
    P00 = P0(2*kk+1)*P0(2*kk+2);
    P01 = P0(2*kk+1)*P1(2*kk+2);
    P10 = P1(2*kk+1)*P0(2*kk+2);
    P11 = P1(2*kk+1)*P1(2*kk+2);
    
    cProb(1, 1, kk+1) = P00;
    cProb(1, 3, kk+1) = P11;
    
    cProb(2, 1, kk+1) = P11;
    cProb(2, 3, kk+1) = P00;
    
    cProb(3, 2, kk+1) = P01;
    cProb(3, 4, kk+1) = P10;
    
    cProb(4, 2, kk+1) = P10;
    cProb(4, 4, kk+1) = P01;
end

% Forward Recursion
aProb(1, 1) = 1;        % Set the initial value of Alpha
for kk = 0:lenBit-1
    aProb(:, kk+2) = transpose(cProb(:, :, kk+1))*aProb(:, kk+1)+1e-300;   % Avoid NaN
    aProb(:, kk+2) = aProb(:, kk+2)/sum(aProb(:, kk+2));
end

% Backward Recursion
bProb(1, end) = 1;    % Set the final value of Beta as [1; 0; 0; 0]
for kk = lenBit:-1:1
    bProb(:, kk) = cProb(:, :, kk)*bProb(:, kk+1)+1e-300;
    bProb(:, kk) = bProb(:, kk)/sum(bProb(:, kk));
end

for kk = 0:lenBit-1
    P_a_0 = transpose(aProb(:, kk+1))*...
        (CnctMtx_a_0.*cProb(:, :, kk+1))*bProb(:, kk+2);
    P_a_1 = transpose(aProb(:, kk+1))*...
        (CnctMtx_a_1.*cProb(:, :, kk+1))*bProb(:, kk+2);
    
    Lapp(kk+1) = log((P_a_0+1e-300)/(P_a_1+1e-300));
    
    cProb_b1(1, 1) = P0(2*kk+1);
    cProb_b1(1, 3) = P1(2*kk+1);
    cProb_b1(2, 1) = P1(2*kk+1);
    cProb_b1(2, 3) = P0(2*kk+1);
    cProb_b1(3, 2) = P0(2*kk+1);
    cProb_b1(3, 4) = P1(2*kk+1);
    cProb_b1(4, 2) = P1(2*kk+1);
    cProb_b1(4, 4) = P0(2*kk+1);
    
    cProb_b2(1, 1) = P0(2*kk+2);
    cProb_b2(1, 3) = P1(2*kk+2);
    cProb_b2(2, 1) = P1(2*kk+2);
    cProb_b2(2, 3) = P0(2*kk+2);
    cProb_b2(3, 2) = P1(2*kk+2);
    cProb_b2(3, 4) = P0(2*kk+2);
    cProb_b2(4, 2) = P0(2*kk+2);
    cProb_b2(4, 4) = P1(2*kk+2);
    
    P_b1_0 = transpose(aProb(:, kk+1))*...
        (CnctMtx_b1_0.*cProb_b1(:, :))*bProb(:, kk+2);
    P_b1_1 = transpose(aProb(:, kk+1))*...
        (CnctMtx_b1_1.*cProb_b1(:, :))*bProb(:, kk+2);
    
    P_b2_0 = transpose(aProb(:, kk+1))*...
        (CnctMtx_b2_0.*cProb_b2(:, :))*bProb(:, kk+2);
    P_b2_1 = transpose(aProb(:, kk+1))*...
        (CnctMtx_b2_1.*cProb_b2(:, :))*bProb(:, kk+2);
    
    Lext21_deintlvr(2*kk+1) = log( (P_b1_0+1e-300) / (P_b1_1+1e-300) );
    Lext21_deintlvr(2*kk+2) = log( (P_b2_0+1e-300) / (P_b2_1+1e-300) );
end
end