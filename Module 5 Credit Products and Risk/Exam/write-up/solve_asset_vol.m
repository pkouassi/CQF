function diff_mse = solve_asset_vol(V0, sigmaV, givenSigmaE)

r = 0.02;
D = 5;
T = 1;

d_1 = (1/(sigmaV*sqrt(T))) * ...
    ( log(V0/D) + (r+0.5*sigmaV^2)*T );
d_2 = d_1 - sigmaV*sqrt(T);

E0 = V0*normcdf(d_1,0,1) - D*exp(-r*T)*normcdf(d_2,0,1);
sigmaE = sigmaV*normcdf(d_1,0,1)*V0/E0;

diff_mse = 10*(E0-3)^3 + (sigmaE-givenSigmaE)^2;

end