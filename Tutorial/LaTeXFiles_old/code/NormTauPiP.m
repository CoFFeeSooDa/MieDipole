
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Inputs:                                                              %
 %   nmax  -> double      : Maximum expansion order                     %
 %   theta -> double      : Polar angle (rad)                           %
 %   order -> string      : Ordering of tables ('normal' or 'reversed') %
 % Outputs:                                                             %
 %   NAng  -> struct array: Normalized Tau, Pi, and P functions         %
 %    .NTau -> double [nx(2n+1)]  : Normalized Tau array                %
 %    .NPi  -> double [nx(2n+1)]  : Normalized Pi array                 %
 %    .NP   -> double [nx(2n+1)]  : Normalized P array                  %
 % Calling functions:                                                   %
 %   Wigner_d                                                           %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NAng = NormTauPiP(nmax,theta,order)
    % Preallocation
    NTau = zeros(nmax,2*nmax+1);
    NPi = zeros(nmax,2*nmax+1);
    NP = zeros(nmax,2*nmax+1);
    for indn = 1:nmax
        % Calling Wigner d Matrix of Order n
        dn = Wigner_d(indn,theta);
        % Setting the Order
        if strcmp(order,'normal') == 1
            dnp1 = dn(:,indn+2)'; % d_(m,+1)^n
            dn01 = dn(:,indn+1)'; % d_(m,0)^n
            dnn1 = dn(:,indn)';   % d_(m,-1)^n
        elseif strcmp(order,'reversed') == 1
            dnp1 = fliplr(dn(:,indn+2)'); % d_(m,+1)^n
            dn01 = fliplr(dn(:,indn+1)'); % d_(m,0)^n
            dnn1 = fliplr(dn(:,indn)');   % d_(m,-1)^n
        end
        % Normalization Constants
        NormTauPi = sqrt((2*indn+1)/8);
        NormP = sqrt((2*indn+1)/2./indn./(indn+1));
        % Output Functions
        NPi(indn,1:2*indn+1) = -NormTauPi.*(dnp1 + dnn1);
        NTau(indn,1:2*indn+1) = -NormTauPi.*(dnp1 - dnn1);
        NP(indn,1:2*indn+1) = NormP.*dn01;
        % Correction to the Floating Numbers
        NPi(abs(NPi)<1e-15) = 0;
        NTau(abs(NTau)<1e-15) = 0;
        NP(abs(NP)<1e-15) = 0;
        % Output a Structure Array
        NAng.NPi = NPi;
        NAng.NTau = NTau;
        NAng.NP = NP;
    end
end