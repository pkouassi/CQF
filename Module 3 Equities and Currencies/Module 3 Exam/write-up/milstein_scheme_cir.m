function r_out = milstein_scheme_cir(r0, eta, gamma, alpha, tau, omega)

    % Define simulation variable
    d_t = tau / 365;
    
    % Initialization of the simulation
    r = zeros(omega, 1+tau/d_t);
    rn = randn(omega, 1+tau/d_t);
    r(:,1) = r0;
    
    % Simulation process
    for sce = 1:omega
        for t = 1:(tau/d_t)
            phi = rn(sce, t);
            r(sce, t+1) = r(sce, t) + (eta - gamma*r(sce, t)) * d_t + ...
                sqrt(alpha*r(sce, t)*d_t)*phi + 0.25*alpha*(phi^2-1)*d_t;
        end
    end
    
    r_out = mean(r,1);
end

