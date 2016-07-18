r = 0.02;
D = 5;
T = 1;
K = 5;

sigmaE = 0.50:0.01:0.80;
Merton_PD = zeros(1, length(sigmaE));
BC_PD = zeros(1, length(sigmaE));

for i = 1:length(sigmaE)
    [results,~] = fsolve(@(x) solve_asset_vol(x(1), x(2), sigmaE(i)),[9;0.25], ...
        optimoptions('fmincon','MaxFunEvals',10000,'MaxIter',10000));
    
    % probability of default Merton
    Merton_PD(i) = normcdf(-((1/(results(2)*sqrt(T))) * ...
        ( log(results(1)/D) + (r+0.5*results(2)^2)*T ) - results(2)*sqrt(T)),0,1);
    
    % probability of default Black-Cox
    h_1 = (log(K/(exp(r*T)*results(1))) + (results(2)^2*T/2))/(results(2)*T);
    h_2 = h_1 - results(2)*sqrt(T);
    BC_PD(i) = normcdf(h_1,0,1) + exp(2*(r-results(2)^2/2) ...
        *log(K/results(2))/results(2)^2)*normcdf(h_2,0,1);
end

avg_diff_PD = mean(BC_PD(sigmaE >= 0.60) - Merton_PD(sigmaE >= 0.60));

plot(sigmaE, Merton_PD, '-b', 'LineWidth', 3);
hold on
plot(sigmaE, BC_PD, '-r', 'LineWidth', 3);
hold off
legend('Merton','Black-Cox')
title('The Comparison on Probability of Default between Merton and Black-Cox')