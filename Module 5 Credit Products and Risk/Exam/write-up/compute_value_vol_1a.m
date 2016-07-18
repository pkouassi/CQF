% solve V0 and sigma_V 
[results,fval] = fsolve(@(x) compute_E0(x(1),x(2)),[9;0.25], ...
    optimoptions('fmincon','MaxFunEvals',10000,'MaxIter',10000));
% check answer
compute_E0(results(1), results(2))
