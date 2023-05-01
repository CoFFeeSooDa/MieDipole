
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Inputs:                                                              %
 %   nr      -> double (1x3): Relative refractive index                 %
 %   ks      -> double      : Dimensionless boundary (k*sphere_radius)  %
 %   nmax    -> double      : Maximum expansion order                   %
 % Outputs:                                                             %
 %   Coeffs  -> struct array: Mie coefficients                          %
 %    .alpha -> double (1xn): Mie coefficient alpha                     %
 %    .beta  -> double (1xn): Mie coefficient beta                      %
 % Calling functions:                                                   %
 %   SphBessel                                                          %
 %   Dlog                                                               %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Coeffs = MieCoreShell(nr,ks,nmax)
    % Checking Input
    if or(ne(max(size(nr)),3) ,ne(max(size(ks)),2) ) == 1
        disp('Error input size of "nr" or "ks" from "MieCoreShell"');
    end
    % Defining Variables
    n0 = nr(1); n1 = nr(2); n2 = nr(3);
    if n0 == n1
        n1 = n1+1e-7;
    end 
    n0kr1=n0*ks(1); n1kr1=n1*ks(1); n1kr2=n1*ks(2); n2kr2=n2*ks(2);
    % Calling Radial Functions (out)
    n0kr1Rad = SphBessel(n0kr1,nmax,1,'bessel') ; n0kr1psi = n0kr1Rad.psi;
    n0kr1Rad = SphBessel(n0kr1,nmax,1,'hankel1'); n0kr1xi  = n0kr1Rad.xi;
    % Calling Radial Functions (shell)
    n1kr1Rad = SphBessel(n1kr1,nmax,1,'hankel1'); n1kr1xi  = n1kr1Rad.xi;
    n1kr2Rad = SphBessel(n1kr2,nmax,1,'hankel1'); n1kr2xi  = n1kr2Rad.xi;
    n1kr1Rad = SphBessel(n1kr1,nmax,1,'bessel') ; n1kr1psi = n1kr1Rad.psi;
    n1kr2Rad = SphBessel(n1kr2,nmax,1,'bessel') ; n1kr2psi = n1kr2Rad.psi;
    % Calling Dlog 
    [n1kr2D1,n1kr2D3]=Dlog(n1kr2,nmax);[n1kr1D1,n1kr1D3]=Dlog(n1kr1,nmax);
    [n0kr1D1,n0kr1D3]=Dlog(n0kr1,nmax);[n2kr2D1,~]=Dlog(n2kr2,nmax);
    % Factors
    f1=n1kr2xi./n1kr2psi; f2=n1kr1xi./n1kr1psi; f3=n0kr1psi./n0kr1xi;
    A  = (n2*n1kr2D3 - n1*n2kr2D1)./(n1*n2kr2D1 - n2*n1kr2D1).*f1;
    B  = (n2*n2kr2D1 - n1*n1kr2D3)./(n1*n1kr2D1 - n2*n2kr2D1).*f1;
    A1 = (n1*n0kr1D1 - n0*n1kr1D3)./(n0*n1kr1D1 - n1*n0kr1D1).*f2;
    A2 = (n0*n1kr1D3 - n1*n0kr1D3)./(n1*n0kr1D3 - n0*n1kr1D1).*f2;
    B1 = (n1*n1kr1D3 - n0*n0kr1D1)./(n0*n0kr1D1 - n1*n1kr1D1).*f2;
    B2 = (n0*n0kr1D3 - n1*n1kr1D3)./(n1*n1kr1D1 - n0*n0kr1D3).*f2;
    % Coefficients
    Coeffs.alpha = (A1-A)./(A2-A).*f3.*...
                   (n0*n1kr1D1 - n1*n0kr1D1)./(n1*n0kr1D3 - n0*n1kr1D1);
    Coeffs.beta  = (B1-B)./(B2-B).*f3.*...
                   (n0*n0kr1D1 - n1*n1kr1D1)./(n1*n1kr1D1 - n0*n0kr1D3);
end