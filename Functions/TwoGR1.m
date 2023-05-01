%% Two-Points Green's Function (in Region 1)

%% Function
function Output = TwoGR1(Settings)
    % Variables
    nmax = Settings.nmax;
    if strcmp(Settings.BC,'simplecavity') == 1
        rhoD = Settings.nr(2)*Settings.k0*Settings.DPos.Sph(1);
    else
        rhoD = Settings.nr(1)*Settings.k0*Settings.DPos.Sph(1);
    end
    rhoA = Settings.nr(2)*Settings.k0*Settings.APos.Sph(1);
    % Radial Functions
    if isfield(Settings,'DRad')  == 0
        if strcmp(Settings.BC,'simplecavity') == 1
            Settings.DRad = SphBessel(rhoD,nmax,1,'bessel');
        else
            Settings.DRad = SphBessel(rhoD,nmax,1,'hankel1');
        end
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
        Settings.Source = SourCoeff(Settings,"Green's function only");
    end
    if isfield(Settings,'layer1') == 0
        if strcmp(Settings.BC,'sphere') == 1
            Settings.Layer1 = MieSingle(Settings.nr,Settings.k0s,nmax);
            Settings.Layer1.d = Settings.Source.p.*transpose(Settings.Layer1.delta);
            Settings.Layer1.c = Settings.Source.q.*transpose(Settings.Layer1.gamma);
        elseif strcmp(Settings.BC,'simplecavity') == 1
            Settings.Layer1 = MieSimCav(Settings.nr,Settings.k0s,nmax);
            Settings.Layer1.d = Settings.Source.r.*transpose(Settings.Layer1.delta);
            Settings.Layer1.c = Settings.Source.s.*transpose(Settings.Layer1.gamma);
        elseif strcmp(Settings.BC,'coreshell') == 1
            % Alert of the Unsupported Function
            Output.error1 = 'EFieldR1 for coreshell structures is not supported yet.';
            Settings.Layer1.gamma = 0;
            Settings.Layer1.delta = 0;
            Settings.Layer1.d = Settings.Source.p.*transpose(Settings.Layer1.delta);
            Settings.Layer1.c = Settings.Source.q.*transpose(Settings.Layer1.gamma);
        end
    end
    % Generating M and N Fields
    Temp.AVSF = VectSphFunc(rhoA,nmax,Settings.ARad,Settings.ANAng,emphi);
    % Donor Dipole Field
    if strcmp(Settings.BC,'simplecavity') == 1
        if isfield(Settings,'EdipS2') == 0
            if isfield(Settings.APos,'Sph2') == 0
                Settings.APos.Sph2 = C2S(Settings.APos.Cart-Settings.DPos.Cart);
            end
            % Field in the Secondary Coordinate
            Settings.EdipS2 = EdipField(Settings.nr(2),Settings.k0,Settings.APos.Sph2,Settings.DOri.Cart);
            % Transforming to the Primary Coordinate
            Settings.EdipS1 = S2S(Settings.EdipS2,(Settings.APos.Sph2(2)-Settings.APos.Sph(2)),0);
        end
    else
        Settings.EdipS1 = [0;0;0];
    end
    % Summing All Order of the Field in Layer 1
    Temp.Layer1M = reshape(sum(Temp.AVSF.M.*Settings.Layer1.c,[1,2]),[3,1]);
    Temp.Layer1N = reshape(sum(Temp.AVSF.N.*Settings.Layer1.d,[1,2]),[3,1]);
    % Two-Points Green's Function (G.Dori, 1/m)
    Output.G = Temp.Layer1M + Temp.Layer1N + Settings.EdipS1...
                /(4*pi*(Settings.nr(2)*Settings.k0)^2);
    % Total Electric Field at the Acceptor Position (dipole moment = 1)
    Output.Etot = 4*pi*(Settings.nr(2)*Settings.k0)^2*Output.G;
    % Total Electric Field at the Acceptor Position (SI)
    epsilon0=8.854187817e-12;
    Output.EtotSI = Settings.Dpstrength*(Settings.nr(2)*Settings.k0)^2....
        /epsilon0*Output.G;
    % Total intensity at the Acceptor position
    Output.Int = norm(Output.Etot).^2;
    % Dipole field
    Output.Edip =  Settings.EdipS1;
    % Etot / Edip
    Output.NEtot = Output.Etot./Settings.EdipS1;
    % Other Additional Outputs
        % Note: Feel Free to Add What You Want!
        % Use the Structure Array to Output Data
        % Example: Output.testAPos = Settings.APos.Sph2;
end
