%% Input for wavelength mode (feel free to change variables)
% Dipole position (micrometer)
    % Donor dipole
    Dx = 0; %Do not change the value! Extremely Important!
    Dy = 0; %Do not change the value! Extremely Important!
    Dz = 0.100;
    % Acceptor dipole
    unit = 0.001; % micrometer
    points = linspace(-120,120,481)*unit;
    [Ax,Az] = meshgrid(points,points);
    Ay = 0;
% Dipole orientation
    % Donor dipole (Cartesian coordinate, [x;y;z])
    Doc = [1;0;0];
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
    % Desired wavelength (microns)
    wavelength = 0.353;
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
        nr(:,2) = 2.0   ;
        % Relative refractive index (Region 2)
        nr(:,3) = sqrt(epsilonmat);
        % Radius of the core and the shell [shell, core] (unit: micron)
        rbc = [0.070, 0.060];
        % Dimensionless radial variable
        k0s = k0*rbc;
    end

%% ------- Do not change the following settings ------- %%
% Creating a mode structure array
    % Mode name 
    Settings.ModeName = 'mapping';
    % Donor position
        % Cartesian coordinate
        Settings.DPos.Cart = [Dx;Dy;Dz];
    % Acceptor position
        % Cartesian coordinate
        AxReshape = reshape(Ax,[1,numel(Ax)]);
        AyReshape = Ay*ones([1,numel(Ax)]);
        AzReshape = reshape(Az,[1,numel(Az)]);
        Settings.APos.Cart = [AxReshape;AyReshape;AzReshape];
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
    % Specific wavelength and wavenumber (data point)
    [~,ii] = min(abs(wavelength-lambda));
    Settings.lambda = lambda(ii);
    Settings.k0 = k0(ii);
    % Dielectric function of the environment
    Settings.nr = nr(ii,:);
    % Radial boundary condition
    Settings.rbc = rbc;
    % Radial dimensionless variable for boundary conditions
    Settings.k0s = k0s(ii,:);

clearvars -except Settings Ar Atheta Aphi Ax Ay Az FilePath ...
    FileName Resize