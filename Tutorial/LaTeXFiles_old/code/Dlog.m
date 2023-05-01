
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Inputs:                                                              %
 %   z     -> double      : Input argument                              %
 %   n     -> double      : Expansion order                             %
 % Outputs:                                                             %
 %   D1    -> double (1xn): Logarithmic derivatives (Ricatti-Bessel)    %
 %   D3    -> double (1xn): Logarithmic derivatives (Ricatti-Hankel)    %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D1,D3] = Dlog(z,n)
    nex = n+15;
    D1 = zeros(1,nex);
    D3 = zeros(1,nex);
    for nn = nex:-1:2
        D1(nn-1) = nn/z - 1/(D1(nn)+nn/z);
    end
    psixi=zeros(nex,1);
    psi0xi0 = (1-exp(2i*z))/2;
    D30 = 1i;
    D10 = cot(z);
    psixi(1) = psi0xi0*(1/z-D10)*(1/z-D30);
    D3(1) = D1(1) + 1i/psixi(1);
    for nn=2:nex
        psixi(nn) = psixi(nn-1)*(nn/z-D1(nn-1))*(nn/z-D3(nn-1));
        D3(nn) = D1(nn) + 1i/psixi(nn);
    end
    D1(n+1:n+15) = [];
    D3(n+1:n+15) = [];
end     