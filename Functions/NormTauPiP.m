%% Special Function: Normalized Vector Spherical Harmonics Array

%% Function
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
            % d_(m,+1)^n
            dnp1 = dn(:,indn+2)';
            % d_(m,0)^n
            dn01 = dn(:,indn+1)';
            % d_(m,-1)^n
            dnn1 = dn(:,indn)';
        elseif strcmp(order,'reversed') == 1
            % d_(m,+1)^n
            dnp1 = fliplr(dn(:,indn+2)');
            % d_(m,0)^n
            dn01 = fliplr(dn(:,indn+1)');
            % d_(m,-1)^n
            dnn1 = fliplr(dn(:,indn)');
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
