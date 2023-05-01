
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Inputs:                                                              %
 % ! Settings -> struct array     : Calculation parameters              %
 % !  .nmax   -> double           : Maximum expansion order             %
 % !  .nr     -> double (1xp)**   : Relative refractive index           %
 % !  .k0     -> double           : modulus of wavevector in vacuum     %
 % !  .k0s    -> double [1x(p-1)] : k0*r_i (r_i = radius of boundary)   %
 % !  .BC     -> string           : Boundary condition                  %
 % !  .DPos   -> double (3x1)     : Donor position                      %
 % !   .Sph   -> double (3x1)     : Presented in a spherical coordinate %
 % !   .Cart  -> double (3x1)     : Presented in a Cartesian coordinate %
 % !  .DOri   -> double (3x1)     : Orientation of the Donor dipole     %
 % !   .Cart  -> double (3x1)     : Presented in a Cartesian coordinate %
 %    .DRad   -> struct array     : Set of radial functions (donor)     %
 %     .h1    -> double (1xn)     : Spherical Bessel functions          %
 %     .raddxi-> double (1xn)     : dxi/kr (See more in SphBessel)      %
 %    .DNAng  -> struct array     : Normalized Tau, Pi, and P (reversed)%
 %     .NTau  -> double [nx(2n+1)]: Normalized Tau array                %
 %     .NPi   -> double [nx(2n+1)]: Normalized Pi array                 %
 %     .NP    -> double [nx(2n+1)]: Normalized P array                  %
 %    .DNAngN -> struct array     : Normalized Tau, Pi, and P (normal)  %
 %     .NTau  -> double [nx(2n+1)]: Normalized Tau array                %
 %     .NPi   -> double [nx(2n+1)]: Normalized Pi array                 %
 %     .NP    -> double [nx(2n+1)]: Normalized P array                  %
 %    .emphi  -> double [nx(2n+1)]: Normalized azimuthal functions      %
 %    .Source -> struct array     : Source Expansion coefficients       %
 %     .p     -> double [nx(2n+1)]: Coefficient p                       %
 %     .q     -> double [nx(2n+1)]: Coefficient q                       %
 %    .Scat   -> struct array     : Mie coefficients                    %
 %     .alpha -> double (1xn)     : Mie coefficient alpha               %
 %     .beta  -> double (1xn)     : Mie coefficient beta                %
 % Outputs:                                                             %
 %   Output   -> struct array     : Storage of output data              %
 %    .Purcell-> double           : Purcell factor at DPos              %
 % Temporary data:                                                      %
 %   Temp     -> struct array     : Storage of temporary data           %
 %    .DVSF   -> struct array     : Vector spherical function (acceptor)%
 %     .M     -> double [nx(2n+1)]: Vector spherical function M         %
 %     .N     -> double [nx(2n+1)]: Vector spherical function N         %
 %    .ScatM  -> double (3x1)     : Scattering contribution (M part)    %
 %    .ScatN  -> double (3x1)     : Scattering contribution (N part)    %
 %                                                                      %
 % !: Variables with "!" is required to make the module work.           %
 % **:                                                                  %
 %  The value of p depends on the spherical scatterer:                  %
 %   Single sphere     -> p = 2                                         %
 %   Core/shell sphere -> p = 3                                         %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Output = PurcellR0(Settings)
    % Variables
    nmax = Settings.nmax;
    rhoD = Settings.nr(1)*Settings.k0*Settings.DPos.Sph(1);
    % Radial Functions for Donor 
    if isfield(Settings,'DRad')  == 0    
        Settings.DRad = SphBessel(rhoD,nmax,1,'hankel1');
    end
    % Angular Functions
    if isfield(Settings,'DNAng') == 0
        Settings.DNAng = NormTauPiP(nmax,Settings.DPos.Sph(2),'reversed');
    end
    if isfield(Settings,'DNAngN') == 0
        Settings.DNAngN = NormTauPiP(nmax,Settings.DPos.Sph(2),'normal');
    end
    % Mie Coefficients
    if isfield(Settings,'Source') == 0
        Settings.Source = SourCoeff(Settings,"Green's function only");
    end
    if isfield(Settings,'Scat') == 0
        if strcmp(Settings.BC,'sphere') == 1
            Settings.Scat = MieSingle(Settings.nr,Settings.k0s,nmax);
        elseif strcmp(Settings.BC,'coreshell') == 1
            Settings.Scat = MieCoreShell(Settings.nr,Settings.k0s,nmax);
        end
        Settings.Scat.a = ...
            Settings.Source.p.*transpose(Settings.Scat.alpha);
        Settings.Scat.b = ...
            Settings.Source.q.*transpose(Settings.Scat.beta);
    end
    % Azimuthal Functions
    if isfield(Settings,'emphi') == 0
        emphi = sqrt(1/2/pi);
    end
    % Generating M and N Fields
    Temp.DVSF=VectSphFunc(rhoD,nmax,Settings.DRad,Settings.DNAngN,emphi);
    % Summing All order of the Scattering Field
    Temp.ScatM = reshape(sum(Temp.DVSF.M.*Settings.Scat.b,[1,2]),[3,1]);
    Temp.ScatN = reshape(sum(Temp.DVSF.N.*Settings.Scat.a,[1,2]),[3,1]);
    % Scattering Part at the Donor Position
    Output.EScat = Temp.ScatM + Temp.ScatN;
    % Purcell Factor
    Output.Purcell = 1 + 6*pi/Settings.k0*...
        imag(transpose(Output.EScat)*Settings.DOri.Sph);
end