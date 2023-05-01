
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
 % !  .APos   -> double (3x1)     : Acceptor position                   %
 % !   .Sph   -> double (3x1)     : Presented in a spherical coordinate %
 % !   .Sph2  -> double (3x1)     : Presented in S2 coordinate          %
 % !   .Cart  -> double (3x1)     : Presented in a Cartesian coordinate %
 % !  .DOri   -> double (3x1)     : Orientation of the Donor dipole     %
 % !   .Cart  -> double (3x1)     : Presented in a Cartesian coordinate %
 %    .DRad   -> struct array     : Set of radial functions (donor)     %
 %     .h1    -> double (1xn)     : Spherical Bessel functions          %
 %     .raddxi-> double (1xn)     : dxi/kr (See more in SphBessel)      %
 %    .ARad   -> struct array     : Set of radial functions (acceptor)  %
 %     .j1    -> double (1xn)     : Spherical Bessel functions          %
 %     .raddpsi> double (1xn)     : dxi/kr (See more in SphBessel)      %
 %    .DNAng  -> struct array     : Normalized Tau, Pi, and P (donor)   %
 %     .NTau  -> double [nx(2n+1)]: Normalized Tau array                %
 %     .NPi   -> double [nx(2n+1)]: Normalized Pi array                 %
 %     .NP    -> double [nx(2n+1)]: Normalized P array                  %
 %    .ANAng  -> struct array     : Normalized Tau, Pi, and P (acceptor)%
 %     .NTau  -> double [nx(2n+1)]: Normalized Tau array                %
 %     .NPi   -> double [nx(2n+1)]: Normalized Pi array                 %
 %     .NP    -> double [nx(2n+1)]: Normalized P array                  %
 %    .emphi  -> double [nx(2n+1)]: Normalized azimuthal functions      %
 %    .Source -> struct array     : Source Expansion coefficients       %
 %     .p     -> double [nx(2n+1)]: Coefficient p                       %
 %     .q     -> double [nx(2n+1)]: Coefficient q                       %
 %    .Scat   -> struct array     : Mie coefficients                    %
 %     .gamma -> double (1xn)     : Mie coefficient gamma               %
 %     .delta -> double (1xn)     : Mie coefficient delta               %
 %    .EdipS1 -> double (3x1)     : EdipField in S1 coordinate          %
 %    .EdipS2 -> double (3x1)     : EdipField in S2 coordinate          %
 % Outputs:                                                             %
 %   Output   -> struct array     : Storage of output data              %
 %    .Etot   -> double (3x1)     : Total electric field at APos        %
 %    .Int    -> double           : Electric field intensity at APos    %
 %    .Edip   -> double (3x1)     : Electric dipole field at APos       %
 %    .NEtot  -> double (3x1)     : Normalized Etot at APos (Etot/Edip) %
 % Temporary data:                                                      %
 %   Temp     -> struct array     : Storage of temporary data           %
 %    .AVSF   -> struct array     : Vector spherical function (acceptor)%
 %     .M     -> double [nx(2n+1)]: Vector spherical function M         %
 %     .N     -> double [nx(2n+1)]: Vector spherical function N         %
 %    .Layer1M-> double (3x1)     : Electric field in layer 1 (M part)  %
 %    .Layer1N-> double (3x1)     : Electric field in layer 1 (N part)  %
 %                                                                      %
 % !: Variables with "!" is required to make the module work.           %
 % **:                                                                  %
 %  The value of p depends on the spherical scatterer:                  %
 %   Single sphere     -> p = 2                                         %
 %   Core/shell sphere -> p = 3                                         %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Output = EFieldR1(Settings)
    % Variables
    nmax = Settings.nmax;
    rhoD = Settings.nr(1)*Settings.k0*Settings.DPos.Sph(1);
    rhoA = Settings.nr(2)*Settings.k0*Settings.APos.Sph(1);
    % Radial Functions
    if isfield(Settings,'DRad')  == 0    
        Settings.DRad = SphBessel(rhoD,nmax,1,'hankel1');
    end
    if isfield(Settings,'ARad') == 0
            Settings.ARad = SphBessel(rhoA,nmax,1,'bessel');
    else
        if isfield(Settings.ARad,'j1') == 0
            Settings.ARad = SphBessel(rhoA,nmax,1,'bessel');
        end
    end
    % Angular Functions
    if isfield(Settings,'DNAng') == 0
        Settings.DNAng = NormTauPiP(nmax,Settings.DPos.Sph(2),'reversed');
    end
    if isfield(Settings,'ANAng') == 0
        Settings.ANAng = NormTauPiP(nmax,Settings.APos.Sph(2),'normal');
    end
    % Azimuthal Functions
    if isfield(Settings,'emphi') == 0
        if Settings.APos.Sph(3) == 0
            % For Speed-Up
            emphi = sqrt(1/2/pi);
        else
            % Setting exp(-inf) = 0 for Useless Array Elements
            m = -inf*ones(nmax,2*nmax+1);
            for ii = 1:nmax
                    m(ii,1:2*ii+1) = -ii:1:ii;
            end
            emphi = sqrt(1/2/pi)*exp(1i*m*Settings.APos.Sph(3));
            % Change exp(-inf) = NaN to Zero
            emphi(isnan(emphi)) = 0;
        end
    end
    % Mie Coefficients
    if isfield(Settings,'Source') == 0
        Settings.Source = SourCoeff(Settings,'dipole');
    end
    if isfield(Settings,'layer1') == 0
        if strcmp(Settings.BC,'sphere') == 1
            Settings.Layer1 = MieSingle(Settings.nr,Settings.k0s,nmax);
        elseif strcmp(Settings.BC,'coreshell') == 1
            % Alert of the Unsupported Function
            Output.error1 = ...
                'EFieldR1 for coreshell structures is not supported yet.';
            Settings.Layer1.gamma = 0;
            Settings.Layer1.delta = 0;
        end
        Settings.Layer1.d = ...
            Settings.Source.p.*transpose(Settings.Layer1.delta);
        Settings.Layer1.c = ...
            Settings.Source.q.*transpose(Settings.Layer1.gamma);
    end
    % Generating M and N Fields
    Temp.AVSF = VectSphFunc(rhoA,nmax,Settings.ARad,Settings.ANAng,emphi);
    % Donor Dipole Field
    if isfield(Settings,'EdipS2') == 0
        if isfield(Settings.APos,'Sph2') == 0
            Settings.APos.Sph2 = ...
                C2S(Settings.APos.Cart-Settings.DPos.Cart);
        end
        % Field in the Secondary Coordinate
        Settings.EdipS2 = EdipField(Settings.nr(1),Settings.k0,...
            Settings.APos.Sph2,Settings.DOri.Cart);
        % Transforming to the Primary Coordinate
        Settings.EdipS1 = S2S(Settings.EdipS2,...
            (Settings.APos.Sph2(2)-Settings.APos.Sph(2)),0);
    end
    % Summing All Order of the Field in Layer 1
    Temp.Layer1M = ...
        reshape(sum(Temp.AVSF.M.*Settings.Layer1.c,[1,2]),[3,1]);
    Temp.Layer1N = ...
        reshape(sum(Temp.AVSF.N.*Settings.Layer1.d,[1,2]),[3,1]);
    % Total Electric Field at the Acceptor Position
    Output.Etot = Temp.Layer1M + Temp.Layer1N;
    % Total intensity at the Acceptor position
    Output.Int = norm(Temp.Layer1M + Temp.Layer1N).^2;
    % Dipole field
    Output.Edip =  Settings.EdipS1;
    % Etot / Edip
    Output.NEtot = Output.Etot./Settings.EdipS1;
    % Other Additional Outputs
        % Note: Feel Free to Add What You Want!
        % Use the Structure Array to Output Data
        % Example: Output.testAPos = Settings.APos.Sph2;
end