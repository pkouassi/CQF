% Define Input Variables
S0 = 100;
tau = 1;
E = 100;
Sigma = 0.2;
r = 0.05;
k = 1;
q = 0;

% Stock Price Generation
S_1000 = milstein_scheme_GBM(S0, tau, Sigma, r, 1000, 0);
S_10000 = milstein_scheme_GBM(S0, tau, Sigma, r, 10000, 0);
S_100000 = milstein_scheme_GBM(S0, tau, Sigma, r, 100000, 0);

price_1000 = asian_option_pricing(S_1000, tau, E, r, k, 1);
price_10000 = asian_option_pricing(S_10000, tau, E, r, k, 1);
price_100000 = asian_option_pricing(S_100000, tau, E, r, k, 1);

% Theoretical option price
% Geometric
q_adj_geo = 0.5*(r+q+Sigma^2/6);
sigma_adj_geo = Sigma/sqrt(3);
d1_geo = (1/(sigma_adj_geo*sqrt(tau)))*(log(S0/E)+...
    (r-q_adj_geo+0.5*sigma_adj_geo^2)*tau);
d2_geo = d1_geo - sigma_adj_geo*sqrt(tau);
asian_call_theo_geo = S0*exp(-q_adj_geo*tau)*normcdf(d1_geo)-...
    E*exp(-r*tau)*normcdf(d2_geo);
asian_put_theo_geo = E*exp(-r*tau)*normcdf(-d2_geo)-...
    S0*exp(-q_adj_geo*tau)*normcdf(-d1_geo);

% Arithmetic
M1 = (exp((r-q)*tau)-1)/((r-q)*tau)*S0;
M2 = 2*exp((2*(r-q)+Sigma^2)*tau)*S0^2/((r-q+Sigma^2)*(2*r-2*q+Sigma^2)*tau^2)+...
    2*S0^2/((r-q)*tau^2)*(1/(2*(r-q)+Sigma^2)-exp((r-q)*tau)/(r-q+Sigma^2));
sigma_adj_arith = sqrt(1/tau * log(M2/M1^2));
d1_arith = (1/(sigma_adj_arith*sqrt(tau)))*(log(M1/E)+...
    (0.5*sigma_adj_arith^2)*tau);
d2_arith = d1_arith - sigma_adj_arith*sqrt(tau);
asian_call_theo_arith = exp(-r*tau)*(M1*normcdf(d1_arith)-E*normcdf(d2_arith));
asian_put_theo_arith = exp(-r*tau)*(E*normcdf(-d2_arith)-M1*normcdf(-d1_arith));