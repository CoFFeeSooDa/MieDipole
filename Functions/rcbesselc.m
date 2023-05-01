%% Riccati-Bessel Functions (Complex Argument)
% Ref: Computation of Special Functions (1996)
%        Authors: Shanjie Zhang, Jianming Jin
%----------------------------------------------
% Called by SphBessel.m

%% Function
function [rcj,rcy,drcj,drcy] = rcbesselc(z,n)
    rcj = zeros(1,n+1);
    rcy = zeros(1,n+1);
    drcj = zeros(1,n+1);
    drcy = zeros(1,n+1);
    NM=n;
    if abs(z) < 1e-60
        disp('ricatti-bessel function precision down');
        drcj(1) = 1; %zeroth order
        rcy = -1.0e300*ones(1,n+1);
        drcy = 1.0e300*ones(1,n+1);
        rcy(1) = -1.0;
        drcy(1) = 0.0;
    else
        rcj(1) = sin(z); %zeroth order
        rcy(1) = -cos(z);
        rcj(2) = rcj(1)/z - cos(z); % first order
        rcy(2) = rcy(1)/z - sin(z);
        rcj0 = rcj(1);
        rcj1 = rcj(2);
        RF0 = rcy(1);
        RF1 = rcy(2);
        for Ky = 3:n+1
            RF2 = (2.0*(Ky-1)-1.0)*RF1/z-RF0;
            if abs(RF2) > 1.0e300
                continue
            end
            rcy(Ky) = RF2;
            RF0 = RF1;
            RF1 = RF2; 
        end
        NMy = Ky-1;
        drcy(1) = sin(z);
        drcy(2) = -rcy(2)/z + rcy(1);
        for Ky = 3:NMy+1
            drcy(Ky) = -(Ky-1)*rcy(Ky)/z + rcy(Ky-1);
        end
        if n >= 2
            M = MSTA1(z,200);
            if M < n
                NM = M;
            else
                M = MSTA2(z,n,15);
            end
            F0 = 0.0;
            F1 = 1.0e-100;
            for K = M+1:-1:1
                F = (2.0*(K-1)+3.0)*F1/z - F0;
                if K <= NM+1
                    rcj (K) = F;
                end
                F0 = F1;
                F1 = F;
            end
            if abs(rcj0) > abs(rcj1)
                CS = rcj0/F;
            else
                CS = rcj1/F0;
            end
            for K = 1:NM+1
                rcj(K) = CS*rcj(K);
            end
        end
        drcj(1) = cos(z);
        drcj(2) = -rcj(2)/z + rcj0;
        for K = 3:NM+1
            drcj(K) = -(K-1)*rcj(K)/z + rcj(K-1);
        end
    end
end
