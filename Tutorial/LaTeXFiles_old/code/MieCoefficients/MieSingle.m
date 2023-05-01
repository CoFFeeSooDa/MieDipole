
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Inputs:                                                              %
 %   nr      -> double (1x2): Relative refractive index                 %
 %   ks      -> double      : Dimensionless boundary (k*sphere_radius)  %
 %   n       -> double      : expansion order                           %
 % Outputs:                                                             %
 %   Coeffs  -> struct array: Mie coefficients                          %
 %    .alpha -> double (1xn): Mie coefficient alpha                     %
 %    .beta  -> double (1xn): Mie coefficient beta                      %
 %    .gamma -> double (1xn): Mie coefficient gamma                     %
 %    .delta -> double (1xn): Mie coefficient delta                     %
 % Calling functions:                                                   %
 %   SphBessel                                                          %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Coeffs = MieSingle(nr,ks,nmax)
    % Checking Input
    if or(ne(max(size(nr)),2) ,ne(max(size(ks)),1) ) == 1
        errmes = 'Error input size of "nr" or "ks" from "MieSingle"';
        dispstat(sprintf(errmes),'keepthis','timestamp');
    end
    % Defining Variables
    n0 = nr(1); n1 = nr(2);
    n0kr1 = n0*ks; n1kr1 = n1*ks;
    % Generating Radial Functions
    n0Rad = SphBessel(n0kr1,nmax,1,'bessel');
    n0psi = n0Rad.psi; n0dpsi = n0Rad.dpsi;
    n0Rad = SphBessel(n0kr1,nmax,1,'hankel1');
    n0xi  = n0Rad.xi ; n0dxi  = n0Rad.dxi ;
    n1Rad = SphBessel(n1kr1,nmax,1,'bessel');
    n1psi = n1Rad.psi; n1dpsi = n1Rad.dpsi;
    % Coefficients
    Coeffs.alpha = -(n1*n0dpsi.*n1psi - n0*n0psi.*n1dpsi)./...
			        (n1*n0dxi.*n1psi - n0*n0xi.*n1dpsi);
    Coeffs.beta  = -(n0*n0dpsi.*n1psi - n1*n0psi.*n1dpsi)./...
		            (n0*n0dxi.*n1psi - n1*n0xi.*n1dpsi);
    Coeffs.gamma = n1*(n0dpsi.*n0xi - n0psi.*n0dxi)./...
			       (n1*n1dpsi.*n0xi - n0*n1psi.*n0dxi);
    Coeffs.delta = n1*(n0dpsi.*n0xi - n0psi.*n0dxi)./...
			       (n0*n1dpsi.*n0xi - n1*n1psi.*n0dxi);
end