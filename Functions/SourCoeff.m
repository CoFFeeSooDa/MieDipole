%% Mie Coefficient for A Dipole Source (cgs unit)

%% Function
function Source = SourCoeff(Settings,type)
    % Variables
    nmax = Settings.nmax;
    % Define the refractive index where the source located to "ni"
    if strcmp(Settings.BC,'simplecavity') == 1
        ni = Settings.nr(2);
    else
        ni = Settings.nr(1);
    end
    kr = ni*Settings.k0*Settings.DPos.Sph(1);
    % Preallocation
    m = -inf*ones(nmax,2*nmax+1);
    % Generate Azimuthal Function
    for ii = 1:nmax
            m(ii,1:2*ii+1) = ii:-1:-ii;
    end
    if Settings.DPos.Sph(3) == 0
        % For Speed-Up
        emphi = sqrt(1/2/pi);
    else
        emphi = sqrt(1/2/pi)*exp(1i*m*Settings.DPos.Sph(3));
        emphi(isnan(emphi)) = 0;
    end
    % Generate N and M Functions
    VSF = VectSphFunc(kr,nmax,Settings.DRad,Settings.DNAng,emphi);
    % Calculate Prefactor
    if strcmp(type,"Green's function only") == 1
        prefactor = 1i*(ni*Settings.k0)*(-1).^m;
    elseif strcmp(type,'dipole') == 1
        % Prefactor from the Green's Function and a Dipole (Gaussian Unit)
        prefactor = 4*pi*1i*(ni*Settings.k0)^3*(-1).^m;
    end
    % Output
    if strcmp(Settings.BC,'simplecavity') == 1
        Source.r = prefactor.*TenCont(VSF.N,Settings.DOri.Sph,[3,1]);
        Source.s = prefactor.*TenCont(VSF.M,Settings.DOri.Sph,[3,1]);
    else
        Source.p = prefactor.*TenCont(VSF.N,Settings.DOri.Sph,[3,1]);
        Source.q = prefactor.*TenCont(VSF.M,Settings.DOri.Sph,[3,1]);
    end
end
