function tran_zeros = zdt(McS)
    [eps, ~] = numden(McS);
    P_num = expand(prod(diag(eps)));
    coefs = flip(coeffs(P_num));
    tran_zeros = roots(coefs);
end