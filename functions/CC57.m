function Code = CC57(Bit)
    lenBit = length(Bit);
    Code = zeros(2*lenBit, 1);
    Bit_zeropad = [0; 0; Bit];
    
    for nn = 1:lenBit
        Code(2*(nn-1)+1) = mod( sum(Bit_zeropad(nn:nn+2)'.*[1 0 1]), 2 );
        Code(2*(nn-1)+2) = mod( sum(Bit_zeropad(nn:nn+2)'.*[1 1 1]), 2 );
    end
end