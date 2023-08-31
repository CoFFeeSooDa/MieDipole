
function epsi = Interpolation(wavelength,wavelengthi,epsii)
    if size(wavelengthi,1) == 1
        epsi = epsii*ones(size(wavelength));
    else
        epsi = interp1(wavelengthi,epsii,wavelength,'pchip');
    end
end