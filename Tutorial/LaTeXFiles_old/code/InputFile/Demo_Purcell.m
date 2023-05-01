%% Input for wavelength mode (feel free to change variables)
% Dipole position (micrometer)
    % Donor dipole
    Dx =  0; %Do not change the value! Extremely Important!
    Dy =  0; %Do not change the value! Extremely Important!
    Dz =  0.010;
% Dipole orientation
    % Donor dipole (Cartesian coordinate, [x;y;z])
    Doc = [1;0;0];
% Largest expansion order (n) in a calculation 
nmax = 30;
% Boundary conditions ('sphere' or 'coreshell')
BC = 'sphere';
    % Calling dielectric data
    excelre = ...
       xlsread('.\DielectricFunction\dielectric function.xlsx','Ag_JPCL');
    epsilonmat = excelre(:,2) + 1i*excelre(:,3);
    % Wavelength (microns) 
    lambda = excelre(:,1)*1e-3;
    %Wavenumber (microns)
    k0 = 2*pi./lambda;
    % Setting conditions
    if strcmp(BC,'sphere') == 1
        % Preallocation
        nr = zeros(size(epsilonmat,1),2);
        % Relative refractive index (Region 0)
        nr(:,1) = 1;
        % Relative refractive index (Region 1)
        nr(:,2) = sqrt(epsilonmat);
        % Radius of the sphere (unit: micron)
        rbc = 0.005;
        % Dimensionless radial variable
        k0s = k0*rbc;
    elseif strcmp(BC,'coreshell') == 1
        % Preallocation
        nr = zeros(size(epsilonmat,1),3);
        % Relative refractive index (Region 0)
        nr(:,1) = 1;
        % Relative refractive index (Region 1)
        nr(:,2) = 2;
        % Relative refractive index (Region 2)
        nr(:,3) = sqrt(epsilonmat);
        % Radius of the core and the shell [shell, core] (unit: micron)
        rbc = [0.070, 0.060];
        % Dimensionless radial variable
        k0s = k0*rbc;
    end

%% Settings for Function "myplot" (Default of Exporting Figures )
fplot.colorstyle = '-k';
fplot.range = [-inf,inf,1e0,1e6];
fplot.yscale = 'log';
%fplot.xlabel = '$\mathrm{Wavenumber}~(\mathrm{cm}^{-1})$';
fplot.xlabel = '$\mathrm{Wavelength}~(\mathrm{nm})$';
fplot.ylabel = '$\mathrm{Purcell~Factor}$';
fplot.subaxis = 1;
fplot.subrange = [-inf,inf,-inf,inf];
fplot.subxlabel = '$\mathrm{Wavelength}~(\mathrm{nm})$';

%% ------- Do not change the following settings ------- %%
% Creating a mode structure array
    % Mode name 
    Settings.ModeName = 'Purcell';
    % Donor position
        % Cartesian coordinate
        Settings.DPos.Cart = [Dx;Dy;Dz];
    % Donor orientation
        % Cartesian coordinate
        Settings.DOri.Cart = Doc;
    % Expansion number
    Settings.nmax = nmax;
    % Type of boundary condition
    Settings.BC = BC;
    % Range of wavelength and wavenumber (data point)
    Settings.lambda = lambda;
    Settings.k0 = k0;
    % Dielectric function of the environment
    Settings.nr = nr;
    % Radial boundary condition
    Settings.rbc = rbc;
    % Radial dimensionless variable for boundary conditions
    Settings.k0s = k0s;

clearvars -except Settings lambda k0 nr k0s FilePath FileName fplot Resize