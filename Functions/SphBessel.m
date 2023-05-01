%% Radial Functions
% Output: 
%       z_1:        Spherical Bessel (Hankel) Functions
%       Z:           Riccati-Bessel (-Hankel) Functions
%       dZ:         Derivative of Z
%       raddZ:     rho*dZ

%% Function
function Rad = SphBessel(kr,nmax,array,type)
    if strcmp(type,'bessel') == 1
        if kr == 0
            z1 = zeros(1,nmax);
            Z = zeros(1,nmax);
            dZ = zeros(1,nmax);
            raddZ = zeros(1,nmax);
            raddZ(1) = 2/3;
            if array == 0
                z1 = z1(nmax);
                Z = Z(nmax);
                dZ = dZ(nmax);
                raddZ = raddZ(nmax);
            end
        else
            [csj,~] = sbesselc(kr,nmax);
            if array == 1
                z1 = csj(2:nmax+1);
%                 z_2 = csj(1:nmax);
%                 lindex = 1:nmax;
            else
                z1 = csj(nmax+1);
%                 z_2 = csj(nmax);
%                 lindex = nmax;
            end
            % Riccati-Bessel Functions and their Derivatives
            [rcj,~,drcj,~] = rcbesselc(kr,nmax);
            if array == 0
                Z = rcj(nmax+1);
                dZ = drcj(nmax+1); 
                raddZ=dZ/kr;
            else
                rcj(1)=[];drcj(1)=[];
                Z = rcj;
                dZ = drcj;
                raddZ=dZ/kr;
            end
        end
        Rad.j1 = z1;    
        Rad.psi = Z;
        Rad.dpsi = dZ;
        Rad.raddpsi = raddZ;
    elseif strcmp(type,'hankel1') == 1
        if kr == 0
            if nmax == 0
                z1 = 1-1i*e300;
                Z = -1i;
                dZ = 1;
            else
                z1 = 0-1i*1e300;
                Z = 0-1i*1e300;
                dZ = 1i*1e300;
                raddZ = 0;
            end
        else
            % Spherical Bessel Functions
            [csj,csy] = sbesselc(kr,nmax);
            if or(numel(csj) < (nmax + 1),numel(csy) < (nmax + 1))
                error('Please decrease the expansion order "n".');
            end
            if array == 1
                z1 = csj(2:nmax+1)+1i*csy(2:nmax+1);
            else
                z1 = csj(nmax+1)+1i*csy(nmax+1);
            end
            % Riccati-Bessel Functions and their Derivatives
            [rcj,rcy,drcj,drcy] = rcbesselc(kr,nmax);
            if array == 0
                Z = rcj(nmax+1) + 1i*rcy(nmax+1);
                dZ = drcj(nmax+1) + 1i*drcy(nmax+1); 
                raddZ=dZ/kr;
            else
                rcj(1)=[];rcy(1)=[];
                drcj(1)=[];drcy(1)=[];
                Z = rcj + 1i*rcy;
                dZ = drcj + 1i*drcy;
                raddZ=dZ/kr;
            end
        end
        Rad.h1 = z1;    
        Rad.xi = Z;
        Rad.dxi = dZ;
        Rad.raddxi = raddZ;
    end
end
