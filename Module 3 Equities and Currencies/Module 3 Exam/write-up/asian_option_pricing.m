function price = asian_option_pricing(S, tau, E, r, k, ifPlotAverage)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    if ifPlotAverage == 1
        n = length(S);

        cum_sum_arith = zeros(size(S));
        cum_sum_geo = zeros(size(S));

        for t = 1:size(S,2)
            cum_sum_arith(:,t) = mean(S(:,1:t),2);
            cum_sum_geo(:,t) = geomean(S(:,1:t),2);
        end
        selected_scenarios = randi([1 n],1,1);
        hold on
        plot(S(selected_scenarios,:))
        plot(cum_sum_arith(selected_scenarios,:),'r')
        plot(cum_sum_geo(selected_scenarios,:),'g')
        title('Stock Price, Arithmetic Average and Geometric Average')
        legend('Stock','Arithmetic','Geometric')
        hold off
    end

    % Asian Option Pricing
    avg_arithmetic = mean(S, 2);
    avg_geometric = geomean(S, 2);

    % fixed strike
    price.call_arithmetic_fixed = mean(exp(-r*tau)*max(avg_arithmetic-E, 0));
    price.call_geometric_fixed = mean(exp(-r*tau)*max(avg_geometric-E, 0));

    price.put_arithmetic_fixed = mean(exp(-r*tau)*max(E-avg_arithmetic, 0));
    price.put_geometric_fixed = mean(exp(-r*tau)*max(E-avg_geometric, 0));

    % floating strike
    price.call_arithmetic_floating = mean(exp(-r*tau)*max(S(:,end)-k*avg_arithmetic, 0));
    price.call_geometric_floating = mean(exp(-r*tau)*max(S(:,end)-k*avg_geometric, 0));

    price.put_arithmetic_floating = mean(exp(-r*tau)*max(k*avg_arithmetic-S(:,end), 0));
    price.put_geometric_floating = mean(exp(-r*tau)*max(k*avg_geometric-S(:,end), 0));

end

