%% Array of Normalized Vector Spherical Functions (M and N)

%% Function
function VSF=VectSphFunc(kr,nmax,Rad,NAng,emphi)
	% Preallocation
    VSF.M = zeros(nmax,2*nmax+1,3);
	VSF.N = zeros(nmax,2*nmax+1,3);
    % Extract Radial Functions
    if isfield(Rad,'h1') == 1
        z1 = Rad.h1;
    elseif isfield(Rad,'j1') == 1
        z1 = Rad.j1;
    end
    if isfield(Rad,'raddxi') == 1
        raddz = Rad.raddxi;
    elseif isfield(Rad,'raddpsi') == 1
        raddz = Rad.raddpsi;
    end
	% Construct the Array of Each Order
	n = (1:nmax)';
    % Construct Radz (j_n(kr)/kr)
    if kr == 0
       Radz = zeros(nmax,1);
       Radz(1) = 1/3;
    else
       Radz = transpose(z1)/kr;
    end
	%M Field
	VSF.M(:,:,2) = 1i.*transpose(z1).*NAng.NPi.*emphi;
	VSF.M(:,:,3) = -transpose(z1).*NAng.NTau.*emphi;
	%N Field
    VSF.N(:,:,1) = Radz.*n.*(n+1).*NAng.NP.*emphi;
	VSF.N(:,:,2) = transpose(raddz).*NAng.NTau.*emphi;
	VSF.N(:,:,3) = 1i.*transpose(raddz).*NAng.NPi.*emphi;
end