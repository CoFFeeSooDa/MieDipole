%% Spherical Bessel Functions (Complex Argument)
% Ref: Computation of Special Functions (1996)
%        Authors: Shanjie Zhang, Jianming Jin
%----------------------------------------------
% Called by SphBessel.m

%% Function
function [csj,csy] = sbesselc(z,n)
    a0 = abs(z);
    nm = n;
    if a0 < 1e-60
        csj = ones([1,n+1])*0;
        csy = ones([1,n+1])*(-1e300);
        csy(1) = 1e0;
        return
    end
    csj = zeros([1,n+1]);
    csj(1) = sin(z)/z;
    csj(2) = (csj(1)-cos(z))/z;
    if n >= 2
        csa = csj(1);
        csb = csj(2);
        m = MSTA1(a0,200);
        if m < n
            nm = m;
        else
            m = MSTA2(a0,n,15);
        end
        cf0 = 0.0;
        cf1 = 1.0-100;
        for k = (m):-1:0
            j = k+1; 
            cf = (2.0*k+3.0)*cf1/z-cf0;
            if k <= nm 
                csj(j)=cf;    
            end
            cf0=cf1;
            cf1=cf;
        end
        if abs(csa) > abs(csb)
            cs = csa/cf;
        else 
            cs = csb/cf0;
        end
        for k = 0:min(nm,n)
            j = k+1;
            csj(j) = cs*csj(j);
        end
    end
    csy=1e200;
    csy(1) = -cos(z)/z;
    csy(2) = (csy(1)-sin(z))/z;
    for k = 2:min(nm,n)
        j = k+1;
        if abs(csj(j-1)) >= abs(csj(j-2))
            csy(j) = (csj(j)*csy(j-1)-1.0/z^2)/csj(j-1);
        else
            csy(j) = (csj(j)*csy(j-2)-(2.0*k-1.0)/z^3)/csj(j-2);
        end
    end
end
