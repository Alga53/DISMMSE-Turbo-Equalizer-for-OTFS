function Code_intlvr = Rect_Interlvr( Code, row_Intlvr, col_intlvr)
intlvrMtx = reshape(Code, row_Intlvr, col_intlvr);
Code_intlvr = reshape(transpose(intlvrMtx), row_Intlvr*col_intlvr, 1);
end