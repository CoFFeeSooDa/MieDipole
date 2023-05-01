%% Coordinate Transformation (Spherical to Spherical)

%% Function
function S1 = S2S(S2,S1t,S1p)
cost = cos(S1t);    sint = sin(S1t);
cosp = cos(S1p);   %sinp = sin(S1p);
T = [  cost    -sint   0;
       sint    cost    0;
       0       0       cosp];
S1 = T*S2;
end
