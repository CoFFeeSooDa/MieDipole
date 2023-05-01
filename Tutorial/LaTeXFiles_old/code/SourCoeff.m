
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Inputs:                                                              %
 %   Settings -> struct array     : Calculation parameters              %
 %    .nmax   -> double           : Maximum expansion order             %
 %    .nr     -> double (1xp)**   : Relative refractive index           %
 %    .k0     -> double           : modulus of wavevector in vacuum     %
 %    .DPos   -> double (3x1)     : Donor position                      %
 %     .Sph   -> double (3x1)     : Presented in a spherical coordinate %
 %   Rad      -> struct array     : Set of radial functions             %
 %    .h1     -> double (1xn)     : Spherical Bessel functions          %
 %    .raddxi -> double (1xn)     : dxi/kr (See more in SphBessel)      %
 %   NAng     -> struct array     : Normalized Tau, Pi, and P funcs.    %
 %    .NTau   -> double [nx(2n+1)]: Normalized Tau array                %
 %    .NPi    -> double [nx(2n+1)]: Normalized Pi array                 %
 %    .NP     -> double [nx(2n+1)]: Normalized P array                  %
 % Outputs:                                                             %
 %   Source   -> struct array     : Source Expansion coefficients       %
 %    .p      -> double [nx(2n+1)]: Coefficient p                       %
 %    .q      -> double [nx(2n+1)]: Coefficient q                       %
 % Calling functions:                                                   %
 %   SphBessel                                                          %
 %   VectSphField                                                       %
 %   TenCont (Tensor Contraction)                                       %
 %                                                                      %
 % **:                                                                  %
 %  The value of p depends on the spherical scatterer:                  %
 %   Single sphere     -> p = 2                                         %
 %   Core/shell sphere -> p = 3                                         %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Source = SourCoeff(Settings,type)
    % Variables
    nmax = Settings.nmax;
    n0 = Settings.nr(1);
    kr = n0*Settings.k0*Settings.DPos.Sph(1);
    % Preallocation
    m = -inf*ones(nmax,2*nmax+1);
    % Generate Azimuthal Function
    for ii = 1:nmax
            m(ii,1:2*ii+1) = ii:-1:-ii;
    end
    if Settings.DPos.Sph(3) == 0
        % For Speed-Up
        emphi = sqrt(1/2/pi);
    else
        emphi = sqrt(1/2/pi)*exp(1i*m*Settings.DPos.Sph(3));
        emphi(isnan(emphi)) = 0;
    end
    % Generate N and M Functions
    VSF = VectSphFunc(kr,nmax,Settings.DRad,Settings.DNAng,emphi);
    % Calculate Prefactor
    if strcmp(type,"Green's function only") == 1
        prefactor = 1i*(n0*Settings.k0)*(-1).^m;
    elseif strcmp(type,'dipole') == 1
        % Prefactor from the Green's Function and a Dipole (Gaussian Unit)
        prefactor = 4*pi*1i*(n0*Settings.k0)^3*(-1).^m;
    end
    % Output
    Source.p = prefactor.*TenCont(VSF.N,Settings.DOri.Sph,[3,1]);
    Source.q = prefactor.*TenCont(VSF.M,Settings.DOri.Sph,[3,1]);
end