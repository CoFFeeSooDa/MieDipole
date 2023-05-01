%% Input for wavelength mode (feel free to change variables)
% Dipole position (micrometer)
    % Donor dipole
    Dx =  0; %Do not change the value! Extremely Important!
    Dy =  0; %Do not change the value! Extremely Important!
    Dz =  0.100;
    % Acceptor dipole
    Ax =  0;
    Ay =  0;
    Az = -0.100;
% Dipole orientation
    % Donor dipole (Cartesian coordinate, [x;y;z])
    Doc = [0;0;1];
    % Acceptor dipole (Cartesian coordinate, [x;y;z])
    Aoc = [0;0;1];
% Largest expansion order (n) in a calculation 
nmax = 70;
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
        rbc = 0.085;
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
fplot.range = [-inf,inf,1e28,1e36];
fplot.yscale = 'log';
fplot.xlabel = '$\mathrm{Wavenumber}~(\mathrm{cm}^{-1})$';
fplot.ylabel = '$\mathrm{Coupling~Factor}~(\mathrm{cm}^{-6})$';
fplot.subaxis = 1;
fplot.subrange = [-inf,inf,-inf,inf];
fplot.subxlabel = '$\mathrm{Wavelength}~(\mathrm{nm})$';
% plot.subylabel = '$\mathrm{Coupling~Factor}~(\mathrm{cm}^{-6})$';

%% ------- Do not change the following settings ------- %%
% Creating a mode structure array
    % Mode name 
    Settings.ModeName = 'wavelength';
    % Donor position
        % Cartesian coordinate
        Settings.DPos.Cart = [Dx;Dy;Dz];
    % Acceptor position
        % Cartesian coordinate
        Settings.APos.Cart = [Ax;Ay;Az];
    % Donor orientation
        % Cartesian coordinate
        Settings.DOri.Cart = Doc;
    % Acceptor orientation
        % Cartesian coordinate
        Settings.AOri.Cart = Aoc;
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