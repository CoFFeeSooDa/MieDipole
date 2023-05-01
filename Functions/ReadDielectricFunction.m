%% Function of reading dielectric function from a .csv file

%  | Wavelength (nm) | Re[Dielectrics] | Im[Dielectrics] |
%  |        .        |        .        |        .        |
%  |        .        |        .        |        .        |
%  |        .        |        .        |        .        |
%  |        .        |        .        |        .        |

function [lambda, epsi] = ReadDielectricFunction(filename)
    file_data = readmatrix(filename);
    lambda = file_data(:,1)*1e-9; % unit: m
    epsi = file_data(:,2) + 1i*file_data(:,3);
end