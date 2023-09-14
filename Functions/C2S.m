%% Cartesian Coordinates to Spherical Coordinates

%% function
function rsph = C2S(rcart)
    % Assigning the Cartesian Components
    x = rcart(1,:); y = rcart(2,:); z = rcart(3,:);
    % Radial Distance
    r = sqrt(x.^2+y.^2+z.^2);
    % Polar Angle
    theta = acos(z./r);
    % Assigning Polar Angle when r = 0
    theta(not(any(r,1))) = 0;
    % Azimuthal Angle
    phi = atan2(y,x);
    % Column Form
    rsph = [r;theta;phi];
end
