%% Logarithmic Derivatives of Riccati-Bessel Functions
% Ref: J. Mod. Opt., 63, 2348-2355 (2016)

%% Function
function [D1,D3] = Dlog(z,n)
    % From experience
    nex = n + floor(abs(1.0478*z + 18.692));
    D1 = zeros(1,nex);
    D3 = zeros(1,n);
    for nn = nex:-1:2
        D1(nn-1) = nn/z - 1/(nn/z + D1(nn));
    end
    D1(1) = (z^2*tan(z)+z-tan(z))/(-z^2+z*tan(z));
    D3(1) = (1i*z^2 - z - 1i)/(z^2 + 1i*z);
    for nn = 2:n
        D3(nn) = -nn/z + 1/(nn/z - D3(nn-1));
    end
    D1(n+1:nex) = [];
end  