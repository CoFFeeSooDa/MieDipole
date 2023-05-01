%% Code of Generalized Mie Theory
% Version: 2.4 (2022.04.29)

% Changes in v2.3
% (1) Fix the error when nmax = 3
% (2) Add a new mode for calculating Purcell factors

% Changes in v2.4
% (1) Redefine variables and their naming

% By Ming-Wei Lee

%% Starting Program
% Clean the Workspace
clear
% File to be Calculated
FilePath = '.\InputFile\'; % Folder Path of Input Files
FileName = 'Demo_MappingMode'; % File Name
% Output Figure Size (value = 0~1)
Resize = 1;

%% Information of initiate a job
tic
dispstat('','init'); % One time only initialization
dispstat(sprintf('Beginning the program...'),'keepthis','timestamp');

%% Loading the Input File
run(append(FilePath,FileName,'.m'));
% Information of the Input File
dispstat(append('The file is imported:'),'keepthis','timestamp');
fprintf(2,append(FileName,'\n'));
% Information of the Using Mode
dispstat(sprintf(append('Using mode: ', Settings.ModeName)),...
    'keepthis','timestamp');
% Information of the Using Structure
dispstat(sprintf(append('Using structure: ',Settings.BC)),...
    'keepthis','timestamp');

%% Pre-Processing (Reducing Computation Time)
% Transforming Coordinate
Settings.DPos.Sph = C2S(Settings.DPos.Cart);
Settings.DOri.Sph = ...
    VecTrans(Settings.DOri.Cart,Settings.DPos.Sph(2:3),'C2S');
% Pre-Calculation of Fixed Variables for Each Mode
if strcmp(Settings.ModeName,'wavelength') == 1
    % Times of the 'for loop'
    Settings.nn = size(Settings.nr,1);
    % Coordinate Transformation
    Settings.APos.Sph = C2S(Settings.APos.Cart);
    Settings.AOri.Sph = ...
        VecTrans(Settings.AOri.Cart,Settings.APos.Sph(2:3),'C2S');
    Settings.APos.Sph2 = C2S(Settings.APos.Cart-Settings.DPos.Cart);
    % Angular Functions
    Settings.DNAng = ...
        NormTauPiP(Settings.nmax,Settings.DPos.Sph(2),'reversed');
    Settings.ANAng = ...
        NormTauPiP(Settings.nmax,Settings.APos.Sph(2),'normal');
elseif strcmp(Settings.ModeName,'angle') == 1
    % Times of the 'for loop'
    Settings.nn = size(Settings.APos.Sph,2);
    % Coordinate Transformation
    Settings.AOri.Sph = ...
        VecTrans(Settings.AOri.Cart,Settings.APos.Sph(2:3),'C2S');
    % Radial Functions
    rhoD = Settings.nr(1)*Settings.k0*Settings.DPos.Sph(1);
    Settings.DRad = SphBessel(rhoD,Settings.nmax,1,'hankel1');
    % Angular Functions
    Settings.DNAng = ...
        NormTauPiP(Settings.nmax,Settings.DPos.Sph(2),'reversed');
    % Source Coefficients
    Settings.Source = SourCoeff(Settings,'dipole');
    if Settings.APos.Sph(1) >= Settings.rbc(1)
        % Scattering Coefficients
        if strcmp(Settings.BC,'sphere') == 1
            Settings.Scat = ...
                MieSingle(Settings.nr,Settings.k0s,Settings.nmax);
        elseif strcmp(Settings.BC,'coreshell') == 1
            Settings.Scat = ...
                MieCoreShell(Settings.nr,Settings.k0s,Settings.nmax);
        end
        Settings.Scat.a = ...
            Settings.Source.p.*transpose(Settings.Scat.alpha);
        Settings.Scat.b = ...
            Settings.Source.q.*transpose(Settings.Scat.beta);
    else 
        % Layer1 Coefficients
        if strcmp(Settings.BC,'sphere') == 1
                Settings.Layer1 = ...
                    MieSingle(Settings.nr,Settings.k0s,Settings.nmax);
            elseif strcmp(Settings.BC,'coreshell') == 1
                str1 = 'Core/shell mapping is not yet fully supported.\n';
                str2 = 'Overwrite E-field in the inner region by zero.\n';
                fprintf(2,str1);
                fprintf(2,str2);
                Settings.Layer1.gamma = 0;
                Settings.Layer1.delta = 0;
        end
        Settings.Layer1.d = ...
            Settings.Source.p.*transpose(Settings.Layer1.delta);
        Settings.Layer1.c = ...
            Settings.Source.q.*transpose(Settings.Layer1.gamma);
    end
    % Radial Functions of the Acceptor
    rhoA = Settings.nr(1)*Settings.k0*Settings.APos.Sph(1);
    Settings.ARad = SphBessel(rhoA,Settings.nmax,1,'hankel1');
elseif strcmp(Settings.ModeName,'mapping') == 1
    % Times of the 'for loop'
    Settings.nn = size(Settings.APos.Cart,2);
    % Coordinate Transformation
    Settings.APos.Sph = C2S(Settings.APos.Cart);
    % Radial Functions
    rhoD = Settings.nr(1)*Settings.k0*Settings.DPos.Sph(1);
    Settings.DRad = SphBessel(rhoD,Settings.nmax,1,'hankel1');
    % Angular Functions
    Settings.DNAng = ...
        NormTauPiP(Settings.nmax,Settings.DPos.Sph(2),'reversed');
    % Source Coefficients
    Settings.Source = SourCoeff(Settings,'dipole');
    % Scattering Coefficients
    if strcmp(Settings.BC,'sphere') == 1
        Settings.Scat = ...
            MieSingle(Settings.nr,Settings.k0s,Settings.nmax);
    elseif strcmp(Settings.BC,'coreshell') == 1
        Settings.Scat = ...
            MieCoreShell(Settings.nr,Settings.k0s,Settings.nmax);
    end
    % Layer1 Coefficients
    if strcmp(Settings.BC,'sphere') == 1
            Settings.Layer1 = ...
                MieSingle(Settings.nr,Settings.k0s,Settings.nmax);
        elseif strcmp(Settings.BC,'coreshell') == 1
            str1 = 'Core/shell mapping is not yet fully supported.\n';
            str2 = 'Overwrite E-field in the inner region by zero.\n';
            fprintf(2,str1);
            fprintf(2,str2);
            Settings.Layer1.gamma = 0;
            Settings.Layer1.delta = 0;
    end
    Settings.Scat.a = ...
        Settings.Source.p.*transpose(Settings.Scat.alpha);
    Settings.Scat.b = ...
        Settings.Source.q.*transpose(Settings.Scat.beta);
    Settings.Layer1.d = ...
        Settings.Source.p.*transpose(Settings.Layer1.delta);
    Settings.Layer1.c = ...
        Settings.Source.q.*transpose(Settings.Layer1.gamma);
elseif strcmp(Settings.ModeName,'Purcell') == 1
    % Times of the 'for loop'
    Settings.nn = size(Settings.nr,1);
    % Angular Functions
    Settings.DNAng = ...
        NormTauPiP(Settings.nmax,Settings.DPos.Sph(2),'reversed');
    Settings.DNAngN = ...
        NormTauPiP(Settings.nmax,Settings.DPos.Sph(2),'normal');
end

%% Preallocation
if strcmp(Settings.ModeName,'Purcell') == 1
    EScat = zeros(Settings.nn,3);
    Purcell = zeros(Settings.nn,1);
else
    Etot = zeros(Settings.nn,3);
    NormEtot = zeros(Settings.nn,3);
    Edip = zeros(Settings.nn,3);
end

%% Main Loop
if strcmp(Settings.ModeName,'wavelength') == 1
    for ii = 1:Settings.nn
        Settings.k0 = k0(ii);
        Settings.nr = nr(ii,:);
        Settings.k0s = k0s(ii,:);
        % Determing which Function is Called by the Acceptor Position 
        if Settings.APos.Sph(1) >= Settings.rbc(1)
            Output = EFieldR0(Settings);
        else
            Output = EFieldR1(Settings);
        end
        Etot(ii,:) = transpose(Output.Etot);
        Edip(ii,:) = transpose(Output.Edip);
        NormEtot(ii,:) = transpose(Output.NEtot);
        % Information
        dispstat(sprintf('Progress: %.2f%%',(ii/Settings.nn)*100),...
            'timestamp');
    end
elseif strcmp(Settings.ModeName,'angle') == 1
    for ii = 1:Settings.nn
        Settings.APos.Cart = [Ax(ii); Ay(ii); Az(ii)];
        Settings.APos.Sph = [Ar; Atheta(ii); Aphi];
        % Determing which Function is Called by the Acceptor Position 
        if Settings.APos.Sph(1) >= Settings.rbc(1)
            Output = EFieldR0(Settings);
        else
            Output = EFieldR1(Settings);
        end
        Etot(ii,:) = transpose(Output.Etot);
        Edip(ii,:) = transpose(Output.Edip);
        % Information
        dispstat(sprintf('Progress: %.2f%%',(ii/Settings.nn)*100),...
            'timestamp');
    end
elseif strcmp(Settings.ModeName,'mapping') == 1
    tmp1 = Settings.APos.Cart;
    tmp2 = Settings.APos.Sph;
    for ii = 1:Settings.nn
        Settings.APos.Cart =tmp1(:,ii);
        Settings.APos.Sph =tmp2(:,ii);
        if Settings.APos.Sph(1) >= Settings.rbc(1)
            Output = EFieldR0(Settings);
        else
            Output = EFieldR1(Settings);
        end
        Etot(ii,:) = transpose(Output.Etot);
        Edip(ii,:) = transpose(Output.Edip);
        % Information
        dispstat(sprintf('Progress: %.2f%%',(ii/Settings.nn)*100),...
            'timestamp');
    end
elseif strcmp(Settings.ModeName,'Purcell') == 1
    for ii = 1:Settings.nn
        Settings.k0 = k0(ii);
        Settings.nr = nr(ii,:);
        Settings.k0s = k0s(ii,:);
        % Calculation of Purcell Factor
        Output = PurcellR0(Settings);
        EScat(ii,:) = transpose(Output.EScat);
        Purcell(ii,:) = transpose(Output.Purcell);
        % Information
        dispstat(sprintf('Progress: %.2f%%',(ii/Settings.nn)*100),...
            'timestamp');
    end
end

%% Output Warnings
if isfield(Output,'error1') ==1
        fprintf(2,append(Output.error1,'\n'));
end

%% Post-Processing
if strcmp(Settings.ModeName,'wavelength') == 1
	% Coupling Factor
    CF = abs(Etot*Settings.AOri.Sph).^2;
    % Coupling Factor along R Direction (Vacuum)
    CFdip = abs(Edip*Settings.AOri.Sph).^2;
    % Setting 0/0 to 0 for Etot/Edip
    NormEtot(isnan(NormEtot)) = 0;
    % Enhancement Factor
    EF = abs(NormEtot*Settings.AOri.Sph).^2;
elseif strcmp(Settings.ModeName,'angle') == 1
    % Coupling Factor along R Direction (Spheres)
    CF = abs(Etot*[1;0;0]).^2;
    % Coupling Factor along R Direction (Vacuum)
    CFdip = abs(Edip*[1;0;0]/sqrt(1)).^2;
elseif strcmp(Settings.ModeName,'mapping') == 1
    % Electric Field Intensity (Spheres)
    EFI = vecnorm(Etot,2,2).^2;
    % Reshape the Array
    EFImap = reshape(EFI,size(Az));
    % Electric Field Intensity (Vacuum)
    EFIdip = vecnorm(Edip,2,2).^2;
    % Reshape the Array
    EFIdipmap = reshape(EFIdip,size(Az));
end

%% Plotting Figures
if strcmp(Settings.ModeName,'wavelength') == 1
    fplot.x = 1./lambda*1e4;
    fplot.y = CF*1e24;
    MyPlot(fplot,Resize,0);
    fplot.y = CFdip*1e24;
    fplot.colorstyle = 'r-';
    MyPlot(fplot,Resize,1);
    fplot.y = EF;
    fplot.colorstyle = 'k-';
    fplot.range = [-inf,inf,1e-2,1e5];
    MyPlot(fplot,Resize,0);
elseif strcmp(Settings.ModeName,'angle') == 1
	fplot.x = Ar*Atheta/Settings.lambda;
    fplot.y = CF*1e24;
    MyPlot(fplot,Resize,0);
    fplot.y = CFdip*1e24;
    fplot.colorstyle = 'b--';
    MyPlot(fplot,Resize,1);
elseif strcmp(Settings.ModeName,'mapping') == 1
    contourf(Ax*1e3,Az*1e3,log10(EFImap*1e24),300,'linestyle','none');
    colormap jet
    colorbar
    hold on
    x = Settings.rbc(1)*linspace(-1,1,101);
    y = sqrt(Settings.rbc(1)^2 - x.^2);
    plot(x*1e3,-y*1e3,'-k','linewidth',2);
    plot(x*1e3,y*1e3,'-k','linewidth',2);
elseif strcmp(Settings.ModeName,'Purcell') == 1
    fplot.x = lambda*1e3;
    fplot.y = Purcell;
    MyPlot(fplot,Resize,0);
end

%% Output information
dispstat('Computation is Finished.','keepprev','timestamp');
toc