% fid = fopen('testjson.json','w');
% fprintf(fid,'%s',jsonencode(k0,PrettyPrint=true));
% 
% 
% fclose(fid);

function result = ReadSettings(filename)
% Initialization
error_msg = false;
% Open the setting file (.json)
fid = fopen(filename);
% Read the file
raw = fread(fid,inf);
% Convert to char
str = char(raw');
% Close the file
fclose(fid);
% Decode
inputjson = jsondecode(str);

% Rename
Settings = inputjson.Settings;
tmp_set = inputjson.tmp_set;
fplot = inputjson.fplot;

% --- post-processing ---
% Default range of the output figure
if isfield(fplot,'range')
    is_nan = isnan(fplot.range);
    for ii = 1:4
        if is_nan(ii) == 1 && (ii == 1 || ii == 3)
            fplot.range(ii) = -inf;
        elseif is_nan(ii) == 1 && (ii == 2 || ii == 4)
            fplot.range(ii) = inf;
        end
    end
end
fplot.range = transpose(fplot.range);
% Default range of the output subfigure
if isfield(fplot,'subrange')
    is_nan = isnan(fplot.subrange);
    for ii = 1:4
        if is_nan(ii) == 1 && (ii == 1 || ii == 3)
            fplot.subrange(ii) = -inf;
        elseif is_nan(ii) == 1 && (ii == 2 || ii == 4)
            fplot.subrange(ii) = inf;
        end
    end
end
fplot.subrange = transpose(fplot.subrange);

% Verify the assignments of dielectric function
if isfield(tmp_set,'epsi0') == 1
    if ischar(tmp_set.epsi0) == 1
        [lambda0, epsi0] = ReadDielectricFunction(tmp_set.epsi0);
    else
        lambda0 = 0; epsi0 = tmp_set.epsi0;
    end
else
    disp("epsi0 isn't assigned yet.");
    error_msg = True;
end
if isfield(tmp_set,'epsi1') == 1
    if ischar(tmp_set.epsi1) == 1
        [lambda1, epsi1] = ReadDielectricFunction(tmp_set.epsi1);
    else
        lambda1 = 0; epsi1 = tmp_set.epsi1;
    end
else
    disp("epsi1 isn't assigned yet.");
    error_msg = True;
end
if strcmp(Settings.BC,'coreshell') == 1
    if isfield(tmp_set,'epsi2') == 1
        if ischar(tmp_set.epsi2) == 1
            [lambda2, epsi2] = ReadDielectricFunction(tmp_set.epsi2);
        else
            lambda2 = 0; epsi2 = tmp_set.epsi2;
        end
    else
        disp("epsi2 isn't assigned yet.");
        error_msg = True;
    end
end

% Post-Processing for Different Mode
if strcmp(Settings.ModeName,'wavelength') == 1
    if strcmp(Settings.BC,'sphere') == 1
        [msize, mind] = max([size(lambda0,1),size(lambda1,1)]);
        nr = zeros(msize,2);
        if mind == 1
            lambda = lambda0; % unit: m
        elseif mind == 2
            lambda = lambda1; % unit: m
        end
        nr(:,1) = sqrt(Interpolation(lambda,lambda0,epsi0));
        nr(:,2) = sqrt(Interpolation(lambda,lambda1,epsi1));
    elseif strcmp(Settings.BC,'coreshell') == 1
        [msize, mind] = max([size(lambda0,1),size(lambda1,1),size(lambda2,1)]);
        nr = zeros(msize,3);
        if mind == 1
            lambda = lambda0; % unit: m
        elseif mind == 2
            lambda = lambda1; % unit: m
        elseif mind == 3
            lambda = lambda2; % unit: m
        end
        nr(:,1) = sqrt(Interpolation(lambda,lambda0,epsi0));
        nr(:,2) = sqrt(Interpolation(lambda,lambda1,epsi1));
        nr(:,3) = sqrt(Interpolation(lambda,lambda2,epsi2));
    end
elseif strcmp(Settings.ModeName,'angle') == 1
    Theta_num = abs(tmp_set.Theta_f - tmp_set.Theta_i)/tmp_set.ThetaResol + 1;
    R = tmp_set.Ar*ones(1,Theta_num);
    Theta = pi*linspace(tmp_set.Theta_i,tmp_set.Theta_f,Theta_num)/180;
    Phi = tmp_set.Phi/180*pi*ones(1,Theta_num);
    Settings.APos.Cart = [tmp_set.Ar*sin(Theta)*cos(tmp_set.Phi);...
                          tmp_set.Ar*sin(Theta)*sin(tmp_set.Phi);...
                          tmp_set.Ar*cos(Theta)];
    Settings.APos.Sph = [R;Theta;Phi];
    if tmp_set.lambda_i == tmp_set.lambda_f
        lambda = tmp_set.lambda_i;
    else
        disp("The initial wavelength isn't equal to the final wavelength.");
    end
    if strcmp(Settings.BC,'sphere') == 1
        nr = zeros(1,2);
        nr(1) = sqrt(Interpolation(lambda,lambda0,epsi0));
        nr(2) = sqrt(Interpolation(lambda,lambda1,epsi1));
    elseif strcmp(Settings.BC,'coreshell') == 1
        nr = zeros(1,3);
        nr(1) = sqrt(Interpolation(lambda,lambda0,epsi0));
        nr(2) = sqrt(Interpolation(lambda,lambda1,epsi1));
        nr(3) = sqrt(Interpolation(lambda,lambda2,epsi2));
    end
elseif strcmp(Settings.ModeName,'mapping') == 1
    Settings.shape = [tmp_set.x_points,tmp_set.y_points];
    % Acceptor dipole
    xgrid = linspace(tmp_set.x_start,tmp_set.x_end,tmp_set.x_points);
    ygrid = linspace(tmp_set.y_start,tmp_set.y_end,tmp_set.y_points);
    if strcmp(tmp_set.plane,'xz') == 1
        [Ax,Az] = meshgrid(xgrid,ygrid);
        Settings.plotx = Ax;
        Settings.ploty = Az;
        
        AxReshape = reshape(Ax,[1,numel(Ax)]);
        AyReshape = tmp_set.third_coord*ones([1,numel(Ax)]);
        AzReshape = reshape(Az,[1,numel(Az)]);
    elseif strcmp(tmp_set.plane,'xy') == 1
        [Ax,Ay] = meshgrid(xgrid,ygrid);
        Settings.plotx = Ax;
        Settings.ploty = Ay;
        
        AxReshape = reshape(Ax,[1,numel(Ax)]);
        AyReshape = reshape(Ay,[1,numel(Ay)]);
        AzReshape = tmp_set.third_coord*ones([1,numel(Ax)]);
    elseif strcmp(tmp_set.plane,'yz') == 1
        [Ay,Az] = meshgrid(xgrid,ygrid);
        Settings.plotx = Ay;
        Settings.ploty = Az;

        AxReshape = tmp_set.third_coord*ones([1,numel(Ay)]);
        AyReshape = reshape(Ay,[1,numel(Ay)]);
        AzReshape = reshape(Az,[1,numel(Az)]);
    end

    Settings.APos.Cart = [AxReshape;AyReshape;AzReshape];
    if tmp_set.lambda_i == tmp_set.lambda_f
        lambda = tmp_set.lambda_i;
    else
        disp("The initial wavelength isn't equal to the final wavelength.");
    end
    if strcmp(Settings.BC,'sphere') == 1
        nr = zeros(1,2);
        nr(1) = sqrt(Interpolation(lambda,lambda0,epsi0));
        nr(2) = sqrt(Interpolation(lambda,lambda1,epsi1));
    elseif strcmp(Settings.BC,'coreshell') == 1
        nr = zeros(1,3);
        nr(1) = sqrt(Interpolation(lambda,lambda0,epsi0));
        nr(2) = sqrt(Interpolation(lambda,lambda1,epsi1));
        nr(3) = sqrt(Interpolation(lambda,lambda2,epsi2));
    end
end


Settings.lambda = lambda;
Settings.k0 = 2*pi./lambda;
Settings.nr = nr;
Settings.rbc = transpose(Settings.rbc);
Settings.k0s = Settings.k0.*Settings.rbc;


result.Settings = Settings;
result.fplot = fplot;
result.error_msg = error_msg;

% Settingstest.lambda = 
end
