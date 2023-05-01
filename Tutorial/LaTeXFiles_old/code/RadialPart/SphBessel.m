
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Inputs:                                                              %
 %   kr        -> double      : Input argument                          %
 %   nmax      -> double      : Maximum expansion order                 %
 %   array     -> double      : Array mode (0 or 1)                     %
 %   type      -> string      : Type of functions ('bessel','hankel1')  %
 % Outputs:                                                             %
 %   Rad       -> struct array: Set of related Bessel (Hankel) functions%
 %    .j1      -> double (1xn): Spherical Bessel functions              %
 %    .psi     -> double (1xn): Riccati-Bessel functions                %
 %    .dpsi    -> double (1xn): Derivatives of psi                      %
 %    .raddpsi -> double (1xn): dpsi/kr                                 %
 %    .h1      -> double (1xn): Spherical Hankel functions              %
 %    .xi      -> double (1xn): Riccati-Hankel functions                %
 %    .dxi     -> double (1xn): Derivatives of xi                       %
 %    .raddxi  -> double (1xn): dxi/kr                                  %
 % Calling functions:                                                   %
 %   sbesselc                                                           %
 %   rcbesselc                                                          %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rad = SphBessel(kr,nmax,array,type)
    if strcmp(type,'bessel') == 1
        if kr == 0
            if nmax == 0
                z1 = 1;
                dZ = 1;
            else
                z1 = 0;
                Z = 0;
                dZ = 0;
                raddZ = 0;
            end
        else
            [csj,~] = sbesselc(kr,nmax);
            if array == 1
                z1 = csj(2:nmax+1);
                z_2 = csj(1:nmax);
                lindex = 1:nmax;
            else
                z1 = csj(nmax+1);
                z_2 = csj(nmax);
                lindex = nmax;
            end
            % Riccati-Bessel Functions
            Z = kr*z1;
            % Derivative of Riccati-Bessel Function
            dZ = kr.*z_2-lindex.*z1;
            raddZ = dZ/kr;
        end
        Rad.j1 = z1; Rad.psi = Z; Rad.dpsi = dZ; Rad.raddpsi = raddZ;
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
        Rad.h1 = z1; Rad.xi = Z; Rad.dxi = dZ; Rad.raddxi = raddZ;
    end
end