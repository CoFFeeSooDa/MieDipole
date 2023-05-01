function result=envj(n,z)
    n = max(1,abs(n));
    result = 0.5*log10(6.28*n)-n*log10(1.36*z/n);
end