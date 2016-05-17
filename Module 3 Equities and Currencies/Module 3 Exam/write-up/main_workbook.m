% Define Input Variables
S0 = 100;
tau = 1;
E = 100;
Sigma = 0.2;
r = 0.05;
k = 1;

% Stock Price Generation
S_1000 = milstein_scheme_GBM(S0, tau, Sigma, r, 1000, 0);
S_10000 = milstein_scheme_GBM(S0, tau, Sigma, r, 10000, 0);
S_100000 = milstein_scheme_GBM(S0, tau, Sigma, r, 100000, 0);

price_1000 = asian_option_pricing(S_1000, tau, E, r, k, 1);
price_10000 = asian_option_pricing(S_10000, tau, E, r, k, 1);
price_100000 = asian_option_pricing(S_100000, tau, E, r, k, 1);
