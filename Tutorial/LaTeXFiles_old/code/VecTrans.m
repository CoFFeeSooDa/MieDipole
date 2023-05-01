function vfin = VecTrans(vini,solidangle,type)
    % Solid Angle to Polar Angle and Azimuthal Angle
    theta = solidangle(1); phi = solidangle(2);
    % Calculating cos(theta) and sin(theta)
    if theta == 0
        sint = 0;
        cost = 1;
    elseif theta == pi
        sint=0;
        cost = -1;
    elseif theta == pi/2 || theta ==3*pi/2
        cost = 0;
        sint = 1;
    else
        cost = cos(theta);
        sint = sin(theta);
    end
    % Calculating cos(phi) and sin(phi)
    if phi == 0
        sinp = 0;
        cosp = 1;
    elseif phi == pi/2
        cosp = 0;
        sinp = 1;
    else
        cosp = cos(phi);
        sinp = sin(phi);
    end
    % Transform Matrix
    T = [ sint*cosp       cost*cosp     -sinp;
            sint*sinp       cost*sinp       cosp;
            cost             -sint              0       ];
    % Asigning the Type of Transformation
    if strcmp(type,'C2S') == 1
        T = transpose(T);
    elseif strcmp(type,'S2C') == 0
        disp('Error from function "VecTrans"');
    end
    vfin = T*vini;
end