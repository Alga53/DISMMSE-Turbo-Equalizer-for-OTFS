function Code_deintlvr = Rect_Deinterlvr( Code, row_Intlvr, col_intlvr)
intlvrMtx = reshape(Code, col_intlvr, row_Intlvr);
Code_deintlvr = reshape(transpose(intlvrMtx), col_intlvr*row_Intlvr, 1);
end