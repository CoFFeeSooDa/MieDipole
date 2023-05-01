%% Two-Points Green's Function (in Region 0)

%% Function
function Output = TwoGR0(Settings)
    % Variables
    nmax = Settings.nmax;
    rhoD = Settings.nr(1)*Settings.k0*Settings.DPos.Sph(1);
    rhoA = Settings.nr(1)*Settings.k0*Settings.APos.Sph(1);
    % Radial Functions
    if isfield(Settings,'DRad')  == 0    
        Settings.DRad = SphBessel(rhoD,nmax,1,'hankel1');
    end
    if isfield(Settings,'ARad')  == 0    
        Settings.ARad = SphBessel(rhoA,nmax,1,'hankel1');
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
            % Speed-Up
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
    % Source Coefficients
    if isfield(Settings,'Source') == 0
        Settings.Source = SourCoeff(Settings,"Green's function only");
    end
    % Mie Coefficients
    if isfield(Settings,'Layer0') == 0
        if strcmp(Settings.BC,'sphere') == 1
            Settings.Layer0 = MieSingle(Settings.nr,Settings.k0s,nmax);
            Settings.Layer0.a = Settings.Source.p.*transpose(Settings.Layer0.alpha);
            Settings.Layer0.b = Settings.Source.q.*transpose(Settings.Layer0.beta);
        elseif strcmp(Settings.BC,'simplecavity') == 1
            Settings.Layer0 = MieSimCav(Settings.nr,Settings.k0s,nmax);
            Settings.Layer0.a = Settings.Source.r.*transpose(Settings.Layer0.alpha);
            Settings.Layer0.b = Settings.Source.s.*transpose(Settings.Layer0.beta);
        elseif strcmp(Settings.BC,'coreshell') == 1
            Settings.Layer0 = MieCoreShell(Settings.nr,Settings.k0s,nmax);
            Settings.Layer0.a = Settings.Source.p.*transpose(Settings.Layer0.alpha);
            Settings.Layer0.b = Settings.Source.q.*transpose(Settings.Layer0.beta);
        end
        Settings.Layer0.a = Settings.Source.p.*transpose(Settings.Layer0.alpha);
        Settings.Layer0.b = Settings.Source.q.*transpose(Settings.Layer0.beta);
    end
    % Generating M and N Fields
    Temp.AVSF=VectSphFunc(rhoA,nmax,Settings.ARad,Settings.ANAng,emphi);
    % Donor Dipole Field
    if strcmp(Settings.BC,'simplecavity') == 1
        Settings.EdipS1 = [0;0;0];
    else
        if isfield(Settings,'EdipS1') == 0
            if isfield(Settings,'EdipS2') == 0
                if isfield(Settings.APos,'Sph2') == 0
                    Settings.APos.Sph2 = C2S(Settings.APos.Cart-Settings.DPos.Cart);
                end
                % Field in the Secondary Coordinate
                Settings.EdipS2 = EdipField(Settings.nr(1),Settings.k0,Settings.APos.Sph2,Settings.DOri.Cart);
            end
            % Transforming to the Primary Coordinate
            Settings.EdipS1 = S2S(Settings.EdipS2,(Settings.APos.Sph2(2)-Settings.APos.Sph(2)),0);
        end
    end
    % Summing All order of the Scattering Field
    Temp.Layer0M = reshape(sum(Temp.AVSF.M.*Settings.Layer0.b,[1,2]),[3,1]);
    Temp.Layer0N = reshape(sum(Temp.AVSF.N.*Settings.Layer0.a,[1,2]),[3,1]);
    % Two-Points Green's Function (G.Dori, 1/m)
    Output.G = Temp.Layer0M + Temp.Layer0N + Settings.EdipS1...
                /(4*pi*(Settings.nr(1)*Settings.k0)^2);
    % Total Electric Field at the Acceptor Position (dipole moment = 1)
    Output.Etot = 4*pi*(Settings.nr(1)*Settings.k0)^2*Output.G;
    % Total Electric Field at the Acceptor Position (SI)
    epsilon0=8.854187817e-12;
    Output.EtotSI = Settings.Dpstrength*(Settings.nr(1)*Settings.k0)^2....
        /epsilon0*Output.G;
    % Total Intensity at the Acceptor Position
    Output.Int = norm(Output.Etot).^2;
    % Dipole Field (dipole moment = 1, Gaussian unit)
    Output.Edip =  Settings.EdipS1;
    % Etot / Edip
    Output.NEtot = Output.Etot./Settings.EdipS1;
    % Other Additional Outputs
        % Note: Feel Free to Add What You Want!
        % Use the Structure Array to Output Data
        % Example: output.testAPos = Settings.APos.Sph2;
end
