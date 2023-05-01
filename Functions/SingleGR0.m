%% Single-Point Green's Function (in Region 0)

%% Function
function Output = SingleGR0(Settings)
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
    if isfield(Settings,'Layer0') == 0
        if strcmp(Settings.BC,'sphere') == 1
            Settings.Layer0 = MieSingle(Settings.nr,Settings.k0s,nmax);
        elseif strcmp(Settings.BC,'coreshell') == 1
            Settings.Layer0 = MieCoreShell(Settings.nr,Settings.k0s,nmax);
        end
        Settings.Layer0.a = Settings.Source.p.*transpose(Settings.Layer0.alpha);
        Settings.Layer0.b = Settings.Source.q.*transpose(Settings.Layer0.beta);
    end
    % Azimuthal Functions
    if isfield(Settings,'emphi') == 0
        emphi = sqrt(1/2/pi);
    end
    % Generating M and N Fields
    Temp.DVSF = VectSphFunc(rhoD,nmax,Settings.DRad,Settings.DNAngN,emphi);
    % Summing All order of the Scattering Field
    Temp.Layer0M = reshape(sum(Temp.DVSF.M.*Settings.Layer0.b,[1,2]),[3,1]);
    Temp.Layer0N = reshape(sum(Temp.DVSF.N.*Settings.Layer0.a,[1,2]),[3,1]);
    % Scattering Part at the Donor Position
    Output.EScat = Temp.Layer0M + Temp.Layer0N;
    % DOri.ImG.Dori
    Output.ImG = Settings.k0/6/pi +...
        imag(transpose(Output.EScat)*Settings.DOri.Sph);
    % Purcell Factor
    Output.Purcell = 6*pi/Settings.k0*Output.ImG;
end
