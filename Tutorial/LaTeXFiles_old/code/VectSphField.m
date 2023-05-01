
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Inputs:                                                              %
 %   kr        -> double             : Dimensionless Radial variable    %
 %   nmax      -> double             : Maximum expansion order          %
 %   Rad       -> struct array**     : Set of radial functions          %
 %    .j1      -> double (1xn)       : Spherical Bessel functions       %
 %    .raddpsi -> double (1xn)       : dpsi/kr (See more in SphBessel)  %
 %    .h1      -> double (1xn)       : Spherical Bessel functions       %
 %    .raddxi  -> double (1xn)       : dxi/kr (See more in SphBessel)   %
 %   NAng      -> struct array       : Normalized Tau, Pi, and P funcs. %
 %    .NTau    -> double [nx(2n+1)]  : Normalized Tau array             %
 %    .NPi     -> double [nx(2n+1)]  : Normalized Pi array              %
 %    .NP      -> double [nx(2n+1)]  : Normalized P array               %
 %   emphi     -> double [nx(2n+1)]  : Array of e^(i*m*phi)             %
 % Outputs:                                                             %
 %   VSF       -> struct array       : Vector spherical functions       %
 %    .M       -> double [nx(2n+1)x3]: Vector spherical function M      %
 %    .N       -> double [nx(2n+1)x3]: Vector spherical function N      %
 %                                                                      %
 % **:                                                                  %
 %  The functions, j1 (raddpsi) and h1 (raddxi), in Rad should not      %
 %  appear simultaneously.                                              %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VSF = VectSphFunc(kr,nmax,Rad,NAng,emphi)
    % Preallocation
    VSF.M = zeros(nmax,2*nmax+1,3);
    VSF.N = zeros(nmax,2*nmax+1,3);
    % Extract Radial Functions
    if isfield(Rad,'h1') == 1
        z1 = Rad.h1;
    elseif isfield(Rad,'j1') == 1
        z1 = Rad.j1;
    end
    if isfield(Rad,'raddxi') == 1
        raddz = Rad.raddxi;
    elseif isfield(Rad,'raddpsi') == 1
        raddz = Rad.raddpsi;
    end
    % Construct the Array of Each Order
    n = (1:nmax)';
    %M Field
    VSF.M(:,:,2) = 1i.*transpose(z1).*NAng.NPi.*emphi;
    VSF.M(:,:,3) = -transpose(z1).*NAng.NTau.*emphi;
    %N Field
    VSF.N(:,:,1) = transpose(z1)/kr.*n.*(n+1).*NAng.NP.*emphi;
    VSF.N(:,:,2) = transpose(raddz).*NAng.NTau.*emphi;
    VSF.N(:,:,3) = 1i.*transpose(raddz).*NAng.NPi.*emphi;
end