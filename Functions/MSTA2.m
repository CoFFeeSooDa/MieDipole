%% Sub-Functions of Spherical Bessel/Riccati-Bessel Functions
% Ref: Computation of Special Functions (1996)
%        Authors: Shanjie Zhang, Jianming Jin
%----------------------------------------------
% Called by sbesselc.m and rcbesselc.m

%% Function
function result = MSTA2(z,n,mp)
    a0 = abs(z);
    hmp = 0.5*mp;
    ejn = envj(n,a0);
    if ejn <= hmp
        obj = mp;
        n0 = fix(1.1*a0);
    else
        obj = hmp+ejn;
        n0 = n;
    end
    f0 = envj(n0,a0)-obj;
    n1 = n0+5;
    f1 = envj(n1,a0)-obj;
    for i = 1:20
        nn = n1-(n1-n0)/(1.0-f0/f1);
        nn = fix(nn); % type conversion: nn should be int
        f = envj(nn,a0)-obj;
        if abs(nn-n1) < 1
           break
        end
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    end
    result = nn+10;
end
