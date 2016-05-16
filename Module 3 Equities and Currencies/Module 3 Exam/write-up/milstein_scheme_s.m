% The stock price simulation under Milstein Scheme
S0 = 100;
tau = 1;
E = 100;
Sigma = 0.2;
r = 0.5;

% Define simulation variable
d_t = 1 / 365;
omega = 1000;

% Initialization of the simulation
S_1 = zeros(omega, 1+tau/d_t);
S_2 = zeros(omega, 1+tau/d_t);
rn = randn(omega, 1+tau/d_t);
S_1(:,1) = S0;
S_2(:,1) = S0;

% Simulation process
for sce = 1:omega
    for t = 1:(tau/d_t)
        phi = rn(sce, t);
        S_1(sce, t+1) = S_1(sce, t) * (1 + r*d_t + Sigma*phi*sqrt(d_t) ...
            + Sigma^2*(phi^2-1)*d_t);
        S_2(sce, t+1) = S_2(sce, t) * (1 + r*d_t + Sigma*(-phi)*sqrt(d_t) ...
            + Sigma^2*(phi^2-1)*d_t);
    end
end

S = [S_1;S_2];