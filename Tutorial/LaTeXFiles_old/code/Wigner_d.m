
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Inputs:                                                              %
 %   j     -> double      : j-dimension SO(3) irreducible representation%
 %   theta -> double      : Polar angle (rad)                           %
 % Outputs:                                                             %
 %   d     -> double (jxj): Wigner d matrix                             %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = Wigner_d(j, theta)
    % Calculation of J+
    m = -j:j-1;
    J = diag(sqrt((j-m).*(j+m+1)), -1);
    % Create the Spectral Decomposition Matrix J_y at z Representation
    Jy = (J-J')/2i;
    % Diagonalization
    [V, D] = eig(Jy);
    % Unitary Transformation
    d = V*diag(exp(-1i*theta*diag(D)))*V';
    % Check the Quality of the Transformation
    if max(max(abs(imag(d)))) > 1e-12
        warn_mes = 'Wigner_d may not give reliable results.';
        dispstat(sprintf(warn_mes),'keepthis','timestamp');
    end
    % Change Data Type (double complex -> double real)
    d = real(d);
end