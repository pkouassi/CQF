% solve alpha 
alpha = fsolve(@(x) correlation_int(x)-0.35,0.85, ...
        optimoptions('fmincon','MaxFunEvals',10000,'MaxIter',10000));
% alpha = 0.8777
correlation_int(0.8777)

% call option price
T = 0.5; K = 120; r = 0; S1 = 90; S2 = 110; sigma1 = 0.3; sigma2 = 0.5;
d2_1 = 1/(sigma1*sqrt(T))*(log(S1/K)+(r-sigma1^2/2)*T);
d2_2 = 1/(sigma2*sqrt(T))*(log(S2/K)+(r-sigma2^2/2)*T);

u1 = 1 - normcdf(-d2_1);
u2 = 1 - normcdf(-d2_2);

B = 1/alpha*log(1+((exp(alpha*u1)-1)*(exp(alpha*u2)-1))/(exp(alpha-1)));