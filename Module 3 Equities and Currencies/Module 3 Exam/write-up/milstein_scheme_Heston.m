function S = milstein_scheme_Heston(S0, tau, Heston, r, omega, ifTestPlot)

    % The stock price simulation under Milstein Scheme
%     S0 = 100;
%     tau = 1;
%     E = 100;
%     Sigma = 0.2;
%     r = 0.05;

    % Define simulation variable
    d_t = tau / 365;
    %omega = 1000;
    
    % Heston model parameters
    kappa = Heston.kappa;
    theta = Heston.theta;
    xi = Heston.xi;
    rho = Heston.rho;
    v0 = Heston.v0;

    % Initialization of the simulation
    S = zeros(omega, 1+tau/d_t);
    v = zeros(omega, 1+tau/d_t);
    
    rn1 = randn(omega, 1+tau/d_t);
    rn2 = randn(omega, 1+tau/d_t);
    
    S(:,1) = S0;
    v(:,1) = v0;

    % Simulation process
    for sce = 1:omega
        for t = 1:(tau/d_t)
            phi1 = rn1(sce, t);
            phi2 = rho * rn1(sce, t) + sqrt(1-rho^2)* rn2(sce, t);
            %phi2 = rn2(sce, t);
            v(sce, t+1) = v(sce, t) + kappa*(theta - v(sce, t))*d_t + ...
                xi*phi2*sqrt(v(sce, t)*d_t) + 0.25*xi^2*(phi2^2-1) * d_t;
            S(sce, t+1) = S(sce, t) * (1 + r*d_t + sqrt(v(sce, t+1))*phi1*sqrt(d_t) ...
                + 0.5*sqrt(v(sce, t+1))^2*(phi1^2-1)*d_t);
        end
    end

    if ifTestPlot == 1
        % Test of the ESG
        selected_scenarios = randi([1 omega],1,1);
        selected_time = randi([1 1+tau/d_t],1,1);
        discount_factors = exp(-r*[0:d_t:tau]);
        S_disc = zeros(size(S));
        for sce = 1:(omega)
            S_disc(sce,:) = S(sce,:).*discount_factors;
        end

        subplot(1,3,1)
        plot(S(selected_scenarios, :))
        xlim([0 tau/d_t])
        title('The Stock Price of a Randomly Selected Path')
        subplot(1,3,2)
        qqplot(rn1(:,selected_time));
        subplot(1,3,3)
        plot(mean(S_disc));
        xlim([0 tau/d_t])
        ylim([99.8 100.2])
        hline = refline([0 100]);
        hline.Color = 'r';
        title('Martingale Test: Discounted Average Stock Price at Each Time Step')
    end
end