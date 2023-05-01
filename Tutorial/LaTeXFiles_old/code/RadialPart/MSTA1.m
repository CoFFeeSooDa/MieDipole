function result = MSTA1(z,mp)
    a0 = abs(z);
    n0 = fix(1.1*a0)+1;  n1 = n0+5;
    f0 = envj(n0,a0)-mp; f1 = envj(n1,a0)-mp;
    for i = 1:20
        nn = n1-(n1-n0)/(1.0-f0/f1);
        nn = fix(nn); % type conversion: nn should be int
        f = envj(nn,a0)-mp;
        if abs(nn-n1) < 1
            break
        end
        n0 = n1; n1 = nn;
        f0 = f1; f1 = f;
    end
    result = nn;
end