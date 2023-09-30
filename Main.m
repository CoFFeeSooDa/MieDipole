%% Code of Generalized Mie Theory
% Version: 2.5 (2023.05.01)

% Changes in v2.3
% (1) Fix the error when nmax = 3
% (2) Add a new mode for calculating Purcell factors

% Changes in v2.4
% (1) Redefine variables and their naming

% Changes in v2.5
% (1) Use .json files to input settings

% By Ming-Wei Lee

%% Starting Program
% Clean the Workspace
clear
% Set the Temperary Path
addpath('./Functions/');												
% File to be Calculated
FilePath = './InputFiles/'; % Folder Path of Input Files
FileName = 'JPCL_11_6796(2020)/Figure5i'; % File Name


% Output Figure Size (value = 0~1)
Resize = 1;

%% Information of Initiate a Job
tic
dispstat('','init'); % One time only initialization
dispstat(sprintf('Beginning the program...'),'keepthis','timestamp');

%% Loading the Input File

% --- json input file ---
Inputfile = ReadSettings(append(FilePath,FileName,'.json'));
Settings = Inputfile.Settings;
fplot = Inputfile.fplot;
if strcmp(Settings.ModeName,'wavelength') == 1
    k0 = Inputfile.Settings.k0;
    k0s = Inputfile.Settings.k0s;
    lambda = Inputfile.Settings.lambda;
    nr = Inputfile.Settings.nr;
elseif strcmp(Settings.ModeName,'angle') == 1
    Ar = Settings.APos.Sph(1,1);
    Atheta = Settings.APos.Sph(2,:);
    Aphi = Settings.APos.Sph(3,1);
    Ax = Settings.APos.Cart(1,:);
    Ay = Settings.APos.Cart(2,:);
    Az = Settings.APos.Cart(3,:);
end
% -----------------------
% Information of the Input File
dispstat(append('The file is imported:'),'keepthis','timestamp');
fprintf(2,append(FileName,'\n'));
% Information of the Using Mode
dispstat(sprintf(append('Using mode: ', Settings.ModeName)),'keepthis','timestamp');
% Information of the Using Structure
dispstat(sprintf(append('Using structure: ',Settings.BC)),'keepthis','timestamp');

%% Checking Whether the File Exists Setting Errors
if strcmp(Settings.BC,'simplecavity') == 1
    if norm(Settings.DPos.Cart) >= Settings.rbc 
        error('Error: The donor dipole should be inside the cavity.');
    end
elseif strcmp(Settings.BC,'sphere') == 1
    if norm(Settings.DPos.Cart) <= Settings.rbc 
        error('Error: The donor dipole should be outside the sphere.');
    end
elseif strcmp(Settings.BC,'coreshell') == 1
    if norm(Settings.DPos.Cart) <= Settings.rbc(1) 
        error('Error: The donor dipole should be outside the shell.');
    end
end

%% Pre-Processing (Reducing Computation Time)
% Transforming Coordinate
Settings.DPos.Sph = C2S(Settings.DPos.Cart);
Settings.DOri.Sph = VecTrans(Settings.DOri.Cart,Settings.DPos.Sph(2:3),'C2S');
% Pre-Calculation of Fixed Variables for Each Mode
if strcmp(Settings.ModeName,'wavelength') == 1
    % Times of the 'for loop'
    Settings.nn = size(Settings.nr,1);
    % Coordinate Transformation
    Settings.APos.Sph = C2S(Settings.APos.Cart);
    Settings.AOri.Sph = VecTrans(Settings.AOri.Cart,Settings.APos.Sph(2:3),'C2S');
    Settings.APos.Sph2 = C2S(Settings.APos.Cart-Settings.DPos.Cart);
    % Angular Functions
    Settings.DNAng = NormTauPiP(Settings.nmax,Settings.DPos.Sph(2),'reversed');
    Settings.ANAng = NormTauPiP(Settings.nmax,Settings.APos.Sph(2),'normal');
elseif strcmp(Settings.ModeName,'angle') == 1
    % Times of the 'for loop'
    Settings.nn = size(Settings.APos.Sph,2);
    % Coordinate Transformation
    Settings.AOri.Sph = VecTrans(Settings.AOri.Cart,Settings.APos.Sph(2:3),'C2S');
    % Radial Functions
    rhoD = Settings.nr(1)*Settings.k0*Settings.DPos.Sph(1);
    Settings.DRad = SphBessel(rhoD,Settings.nmax,1,'hankel1');
    % Angular Functions
    Settings.DNAng = NormTauPiP(Settings.nmax,Settings.DPos.Sph(2),'reversed');
    % Source Coefficients
    Settings.Source = SourCoeff(Settings,"Green's function only");
    if Settings.APos.Sph(1) >= Settings.rbc(1)
        % Layer0 Coefficients
        if strcmp(Settings.BC,'sphere') == 1
            Settings.Layer0 = MieSingle(Settings.nr,Settings.k0s,Settings.nmax);
        elseif strcmp(Settings.BC,'coreshell') == 1
            Settings.Layer0 = MieCoreShell(Settings.nr,Settings.k0s,Settings.nmax);
        end
        Settings.Layer0.a = Settings.Source.p.*transpose(Settings.Layer0.alpha);
        Settings.Layer0.b = Settings.Source.q.*transpose(Settings.Layer0.beta);
    else 
        % Layer1 Coefficients
        if strcmp(Settings.BC,'sphere') == 1
                Settings.Layer1 = MieSingle(Settings.nr,Settings.k0s,Settings.nmax);
            elseif strcmp(Settings.BC,'coreshell') == 1
                fprintf(2,'The feature of core/shell mapping is not supported yet.\n');
                fprintf(2,'Overwrite the electric field of the inner region by zero.\n');
                Settings.Layer1.gamma = 0;
                Settings.Layer1.delta = 0;
        end
        Settings.Layer1.d = Settings.Source.p.*transpose(Settings.Layer1.delta);
        Settings.Layer1.c = Settings.Source.q.*transpose(Settings.Layer1.gamma);
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
    if strcmp(Settings.BC,'simplecavity') == 1
        Settings.DRad = SphBessel(rhoD,Settings.nmax,1,'bessel');
    else
        Settings.DRad = SphBessel(rhoD,Settings.nmax,1,'hankel1');
    end
    % Angular Functions
    Settings.DNAng = NormTauPiP(Settings.nmax,Settings.DPos.Sph(2),'reversed');
    % Source Coefficients
    Settings.Source = SourCoeff(Settings,"Green's function only");
    % Layer0 Coefficients
    if strcmp(Settings.BC,'sphere') == 1
        Settings.Layer0 = MieSingle(Settings.nr,Settings.k0s,Settings.nmax);
    elseif strcmp(Settings.BC,'coreshell') == 1
        Settings.Layer0 = MieCoreShell(Settings.nr,Settings.k0s,Settings.nmax);
    elseif strcmp(Settings.BC,'simplecavity') == 1
        Settings.Layer0 = MieSimCav(Settings.nr,Settings.k0s,Settings.nmax);
    end
    % Layer1 Coefficients
    if strcmp(Settings.BC,'sphere') == 1
        Settings.Layer1 = MieSingle(Settings.nr,Settings.k0s,Settings.nmax);
    elseif strcmp(Settings.BC,'coreshell') == 1
        fprintf(2,'The feature of core/shell mapping is not supported yet.\n');
        fprintf(2,'Overwrite the electric field of the inner region by zero.\n');
        Settings.Layer1.gamma = 0;
        Settings.Layer1.delta = 0;
    elseif strcmp(Settings.BC,'simplecavity') == 1
        Settings.Layer1 = MieSimCav(Settings.nr,Settings.k0s,Settings.nmax);
    end
    if strcmp(Settings.BC,'simplecavity') == 1
        Settings.Layer0.a = Settings.Source.r.*transpose(Settings.Layer0.alpha);
        Settings.Layer0.b = Settings.Source.s.*transpose(Settings.Layer0.beta);
        Settings.Layer1.d = Settings.Source.r.*transpose(Settings.Layer1.delta);
        Settings.Layer1.c = Settings.Source.s.*transpose(Settings.Layer1.gamma);
    else
        Settings.Layer0.a = Settings.Source.p.*transpose(Settings.Layer0.alpha);
        Settings.Layer0.b = Settings.Source.q.*transpose(Settings.Layer0.beta);
        Settings.Layer1.d = Settings.Source.p.*transpose(Settings.Layer1.delta);
        Settings.Layer1.c = Settings.Source.q.*transpose(Settings.Layer1.gamma);
    end
end

%% Preallocation
if strcmp(Settings.ModeName,'wavelength') == 1
    EScat = zeros(Settings.nn,3);
    ImG = zeros(Settings.nn,1);
    if Settings.APos.Cart ~= Settings.DPos.Cart
        ImG_vec = zeros(Settings.nn,3);
    end
    Purcell = zeros(Settings.nn,1);
else
    Etot = zeros(Settings.nn,3);
    NormEtot = zeros(Settings.nn,3);
    Edip = zeros(Settings.nn,3);
    EtotSI = zeros(Settings.nn,3);
end

%% Main Loop
if strcmp(Settings.ModeName,'wavelength') == 1
    for ii = 1:Settings.nn
        Settings.k0 = k0(ii);
        Settings.nr = nr(ii,:);
        Settings.k0s = k0s(ii,:);
        % Determing which Function is Called by the Acceptor Position
        if Settings.APos.Cart == Settings.DPos.Cart
            if strcmp(Settings.BC,'simplecavity') == 1
                Output = SingleGR1(Settings);
            else
                Output = SingleGR0(Settings);
            end
            EScat(ii,:) = transpose(Output.EScat);
            ImG(ii,:) = transpose(Output.ImG);
            Purcell(ii,:) = transpose(Output.Purcell);
        else
            if Settings.APos.Sph(1) >= Settings.rbc(1)
                Output = TwoGR0(Settings);
            else
                Output = TwoGR1(Settings);
            end
            ImG_vec(ii,:) = transpose(imag(Output.G));
            Etot(ii,:) = transpose(Output.Etot);
            Edip(ii,:) = transpose(Output.Edip);
            NormEtot(ii,:) = transpose(Output.NEtot);
        end
        % Information
        dispstat(sprintf('Progress: %.2f%%',(ii/Settings.nn)*100),'timestamp');
    end
elseif strcmp(Settings.ModeName,'angle') == 1
    for ii = 1:Settings.nn
        Settings.APos.Cart = [Ax(ii); Ay(ii); Az(ii)];
        Settings.APos.Sph = [Ar; Atheta(ii); Aphi];
        % Determing which Function is Called by the Acceptor Position 
        if Settings.APos.Sph(1) >= Settings.rbc(1)
            Output = TwoGR0(Settings);
        else
            Output = TwoGR1(Settings);
        end
        Etot(ii,:) = transpose(Output.Etot);
        Edip(ii,:) = transpose(Output.Edip);
        % Information
        dispstat(sprintf('Progress: %.2f%%',(ii/Settings.nn)*100),'timestamp');
    end
elseif strcmp(Settings.ModeName,'mapping') == 1
    tmp1 = Settings.APos.Cart;
    tmp2 = Settings.APos.Sph;
    for ii = 1:Settings.nn
        Settings.APos.Cart =tmp1(:,ii);
        Settings.APos.Sph =tmp2(:,ii);
        if Settings.APos.Sph(1) >= Settings.rbc(1)
            Output = TwoGR0(Settings);
        else
            Output = TwoGR1(Settings);
        end
        Etot(ii,:) = transpose(Output.Etot);
        EtotSI(ii,:) = transpose(Output.EtotSI);
        Edip(ii,:) = transpose(Output.Edip);
        % Information
        dispstat(sprintf('Progress: %.2f%%',(ii/Settings.nn)*100),'timestamp');
    end
end

%% Output Warnings
if isfield(Output,'error1') ==1
        fprintf(2,append(Output.error1,'\n'));
end

%% Post-Processing
if strcmp(Settings.ModeName,'wavelength') == 1
    if strcmp(Settings.Quantity,'CF') == 1
        % Coupling Factor
        CF = abs(Etot*Settings.AOri.Sph).^2;
        % Coupling Factor along R Direction (Vacuum)
        CFdip = abs(Edip*Settings.AOri.Sph).^2;
        % Setting 0/0 to 0 for Etot/Edip
        NormEtot(isnan(NormEtot)) = 0;
        % Enhancement Factor
        EF = abs(NormEtot*Settings.AOri.Sph).^2;
    elseif strcmp(Settings.Quantity,'Purcell') == 1

    elseif strcmp(Settings.Quantity,'ImG') == 1
        if exist('ImG_vec','var') == 1
            ImG = ImG_vec*Settings.AOri.Sph;
        end
    elseif strcmp(Settings.Quantity,'J') == 1
        if exist('ImG_vec','var') == 1
            ImG = ImG_vec*Settings.AOri.Sph;
        end
        c = 2.9979e8;
        Debye = 3.33564e-30;
        epsilon0 = 8.854187817e-12;
        hbar = 1.05457182e-34;
        const = (2*pi*1239.84193./(lambda*1e9)*2.4179893e14).^2....
            /c^2*Debye^2/(pi*hbar*epsilon0);
        J = const.*ImG;
    end
elseif strcmp(Settings.ModeName,'angle') == 1
    % Coupling Factor along R Direction (Spheres)
    CF = abs(Etot*[1;0;0]).^2;
    % Coupling Factor along R Direction (Vacuum)
    CFdip = abs(Edip*[1;0;0]/sqrt(1)).^2;
elseif strcmp(Settings.ModeName,'mapping') == 1
%     c=2.9979e8;
%     const = Dpstrength*c^2*1e-5;
    % Electric Field Intensity (Spheres)
    EFI = vecnorm(EtotSI,2,2).^2;
    % Reshape the Array
    EFImap = reshape(EFI,Settings.shape);
    % Electric Field Intensity (Vacuum)
    EFIdip = vecnorm(Edip,2,2).^2;
    % Reshape the Array
    EFIdipmap = reshape(EFIdip,Settings.shape);
end

%% Plotting Figures
if strcmp(Settings.ModeName,'wavelength') == 1
    if strcmp(Settings.Quantity,'CF') == 1
        fplot.x = 1./lambda*1e-2;
        fplot.y = CF*1e-12;
        MyPlot(fplot,Resize,0);
        fplot.y = CFdip*1e-12;
        fplot.colorstyle = 'r-';
        MyPlot(fplot,Resize,1);
        if strcmp(Settings.BC,'sphere') == 1
            legend({'Single Sphere','Vacuum (QED)'},'interpreter','latex');
        elseif strcmp(Settings.BC,'coreshell') == 1
            legend({'Core/Shell Sphere','Vacuum (QED)'},'interpreter','latex');
        end
        fplot.y = EF;
        fplot.colorstyle = '-';
        fplot.range = [-inf,inf,1e-3,1e5];
        fplot.ylabel = 'Enhancement';
        MyPlot(fplot,Resize,0);
        if strcmp(Settings.BC,'sphere') == 1
            legend({'Single Sphere'},'interpreter','latex');
        elseif strcmp(Settings.BC,'coreshell') == 1
            legend({'Core/Shell Sphere'},'interpreter','latex');
        end
    elseif strcmp(Settings.Quantity,'Purcell') == 1
%         fplot.x = 1./lambda*1e4;
        fplot.x = 1239.84193./(lambda*1e9);
        fplot.y = Purcell;    
        MyPlot(fplot,Resize,0);
    elseif strcmp(Settings.Quantity,'ImG') == 1
        fplot.x = 1239.84193./(lambda*1e9);
        fplot.y = ImG;    
        MyPlot(fplot,Resize,0);
    elseif strcmp(Settings.Quantity,'J') == 1
        fplot.x = 1239.84193./(lambda*1e9);
        fplot.y = J;    
        MyPlot(fplot,Resize,0);
    end
elseif strcmp(Settings.ModeName,'angle') == 1
    if strcmp(Settings.Quantity,'CF') == 1
        fplot.x = Ar*Atheta/Settings.lambda;
        fplot.y = CF*1e-12;
        MyPlot(fplot,Resize,0);
        fplot.y = CFdip*1e-12;
        fplot.colorstyle = 'b--';
        MyPlot(fplot,Resize,1);
    end
elseif strcmp(Settings.ModeName,'mapping') == 1
    if strcmp(Settings.Quantity,'CF') == 1
        figure
        contourf(Settings.plotx*1e9,Settings.ploty*1e9,log10(EFImap),300,'linestyle','none');
        colormap jet
        colorbar
        hold on
        x = Settings.rbc(1)*linspace(-1,1,201);
        y = sqrt(Settings.rbc(1)^2 - x.^2);
        plot(x*1e9,-y*1e9,'-k','linewidth',2);
        plot(x*1e9,y*1e9,'-k','linewidth',2);
        figure
        contourf(Settings.plotx*1e9,Settings.ploty*1e9,log10(EFIdipmap),300,'linestyle','none');
        colormap jet
        colorbar
        hold on
        x = Settings.rbc(1)*linspace(-1,1,201);
        y = sqrt(Settings.rbc(1)^2 - x.^2);
        plot(x*1e9,-y*1e9,'-k','linewidth',2);
        plot(x*1e9,y*1e9,'-k','linewidth',2);
    end
    
end

%% Output Information
dispstat('Computation is Finished.','keepprev','timestamp');
% rmpath('./Functions/');
toc
