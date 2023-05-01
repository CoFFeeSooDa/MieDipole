%% Sub-Functions of MSTA1/MSTA2
% Ref: Computation of Special Functions (1996)
%        Authors: Shanjie Zhang, Jianming Jin
%----------------------------------------------
% Called by MSTA1.m and MSTA2.m

%% Function
function result=envj(n,z)
    n = max(1,abs(n));
    result = 0.5*log10(6.28*n)-n*log10(1.36*z/n);
end
