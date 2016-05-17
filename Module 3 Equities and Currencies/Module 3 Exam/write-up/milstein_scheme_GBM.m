function S = milstein_scheme_GBM(S0, tau, Sigma, r, omega, ifTestPlot)

    % The stock price simulation under Milstein Scheme
%     S0 = 100;
%     tau = 1;
%     E = 100;
%     Sigma = 0.2;
%     r = 0.05;

    % Define simulation variable
    d_t = tau / 365;
    %omega = 1000;

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

    if ifTestPlot == 1
    % Test of the ESG
    selected_scenarios = randi([1 2*omega],1,1);
    selected_time = randi([1 1+tau/d_t],1,1);
    discount_factors = exp(-r*[0:d_t:tau]);
    S_disc = zeros(size(S));
    for sce = 1:(2*omega)
        S_disc(sce,:) = S(sce,:).*discount_factors;
    end

    subplot(1,3,1)
        plot(S_1(selected_scenarios, :))
        xlim([0 tau/d_t])
        title('The Stock Price of a Randomly Selected Path')
        subplot(1,3,2)
        qqplot(rn(:,selected_time));
        subplot(1,3,3)
        plot(mean(S_disc));
        xlim([0 tau/d_t])
        ylim([99.8 100.2])
        hline = refline([0 100]);
        hline.Color = 'r';
        title('Martingale Test: Discounted Average Stock Price at Each Time Step')
    end
end