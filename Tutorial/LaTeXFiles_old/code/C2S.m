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
    phi = atan(y./x);
    % Assigning Azimuthal Angle when x = 0 and y = 0
    phi(and(not(any(y,1)),not(any(x,1)))) = 0;
    % Assigning Azimuthal Angle when x < 0
    phi(and(x<0,any(y,1))) = pi;
    % Column Form
    rsph = [r;theta;phi];
end