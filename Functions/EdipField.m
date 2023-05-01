%% Dipole Field under the Secondary Spherical Coordinate (cgs unit)

%% Function
function EdipS = EdipField(nr,k,rdip,vdip)
% Preallocation
NX = zeros(3,1); NY = zeros(3,1); NZ = zeros(3,1);
r = rdip(1); theta = rdip(2); phi = rdip(3);
% Radial Function (1i*k^3 is canceled)
Rad1 = exp(1i*k*r)/r*(r^(-2)-1i*k/r);
Rad2 = exp(1i*k*r)/r*(k^2+1i*k/r-r^(-2));
% Z-Component (1i*k^3 is canceled)
NZ(1) =  Rad1*cos(theta)*2;
NZ(2) = -Rad2*sin(theta);
% X-Component (1i*k^3 is canceled)
NX(1) =  Rad1*sin(theta)*cos(phi)*2;
NX(2) =  Rad2*cos(theta)*cos(phi);
NX(3) = -Rad2*sin(phi);
% Y-Component (1i*k^3 is canceled)
NY(1) =  Rad1*sin(theta)*sin(phi)*2;
NY(2) =  Rad2*cos(theta)*sin(phi);
NY(3) =  Rad2*cos(phi);
% Electric Dipole Field (Gaussian Unit)
EdipS = (NX*vdip(1) + NY*vdip(2) + NZ*vdip(3))*nr;
% EdipC = VecTrans(EdipS,[theta,phi],'S2C');
end
