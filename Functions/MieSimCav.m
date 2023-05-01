%% Mie coefficients for A Simple Spherical Cavity

%% Function
function Coeffs = MieSimCav(nr,ks,nmax)
    % Checking Input
    if or(ne(max(size(nr)),2) ,ne(max(size(ks)),1) ) == 1
        errmes = 'Error input size of "nr" or "ks" from "MieSingle"';
        dispstat(sprintf(errmes),'keepthis','timestamp');
    end
    % Defining Variables
    n0 = nr(1); n1 = nr(2);
    n0kr1 = n0*ks; n1kr1 = n1*ks;
    % Generating Radial Functions
	n0Rad = SphBessel(n0kr1,nmax,1,'hankel1');
    n0xi = n0Rad.xi;
	n1Rad = SphBessel(n1kr1,nmax,1,'bessel');
    n1psi = n1Rad.psi;
    n1Rad = SphBessel(n1kr1,nmax,1,'hankel1');
    n1xi = n1Rad.xi;
    % Calling Dlog
    [n1kr1D1,n1kr1D3] = Dlog(n1kr1,nmax);
    [~,n0kr1D3] = Dlog(n0kr1,nmax);
    % Coefficients
    Coeffs.alpha = (n0*n1kr1D3 - n0*n1kr1D1)./...
			(n1*n0kr1D3 - n0*n1kr1D1).*n1xi./n0xi;
	Coeffs.beta = (n0*n1kr1D3 - n0*n1kr1D1)./...
		    (n0*n0kr1D3 - n1*n1kr1D1).*n1xi./n0xi;
	Coeffs.gamma =  - (n1*n1kr1D3 - n0*n0kr1D3)./...
			(n1*n1kr1D1 - n0*n0kr1D3).*n1xi./n1psi;
	Coeffs.delta =  - (n0*n1kr1D3 - n1*n0kr1D3)./...
			(n0*n1kr1D1 - n1*n0kr1D3).*n1xi./n1psi;
end
