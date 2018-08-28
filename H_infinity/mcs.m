function McS = mcs(P)
    [~, den] = numden(P);
    Den = gcd(den);
    Pp = P*Den;
    S = smithForm(Pp);
    McS = simplify(S/Den);
end