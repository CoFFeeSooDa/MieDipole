 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Inputs:                                                              %
 %   S2    -> double (3x1): Secondary spherical coordinate              %
 %   S1t   -> double      : Theta in the S1 coordinate                  %
 %   S1p   -> double      : Phi in the S1 coordinate                    %
 % Outputs:                                                             %
 %   S1    -> double (3x1): Primary spherical coordinate                %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S1 = S2S(S2,S1t,S1p)
    cost = cos(S1t);
    sint = sin(S1t);
    cosp = cos(S1p);
    % Coordinate Transformation Matrix [Eq.(87)]
    R = [cost    -sint   0;
         sint     cost   0;
         0        0      cosp];
    S1 = R*S2;
end