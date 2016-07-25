% main function for the final projection

%% Part II Credit Value Adjustment
% Create RateSpec from the Interest Rate Curve
% define workplace
crt_path = which('main_cva.m');
parrent_path = strrep(crt_path, 'Code\main_cva.m', '');
data_file_path = strcat(parrent_path, 'Data\');

% read data from directory
swap_rate_raw = xlsread(strcat(data_file_path, 'us_swap_rate.xlsx'));
swap_tenor = swap_rate_raw(:,1);
swap_rates = swap_rate_raw(:,2);

% Read swaps from spreadsheet
swapFile = 'cva-swap-portfolio.xls';
swapData = readtable(strcat(data_file_path,swapFile),'Sheet','Swap Portfolio');
%swapData = readtable(swapFile,'Sheet','Swap Portfolio');

swaps = struct( ...
    'Counterparty',[], ...
    'NettingID',[], ...
    'Principal',[], ...
    'Maturity',[], ...
    'LegRate',[], ...
    'LegType',[], ...
    'LatestFloatingRate',[], ...
    'FloatingResetDates',[]);

swaps.Counterparty = swapData.CounterpartyID;
swaps.NettingID = swapData.NettingID;
swaps.Principal = swapData.Principal;
swaps.Maturity = swapData.Maturity;
swaps.LegType = [swapData.LegType ~swapData.LegType];
swaps.LegRate = [swapData.LegRateReceiving swapData.LegRatePaying];
swaps.LatestFloatingRate = swapData.LatestFloatingRate;
swaps.Period = swapData.Period;
swaps.LegReset = ones(size(swaps.LegType));

numSwaps = numel(swaps.Counterparty);

% Create RateSpec from the Interest Rate Curve
% define workplace
% crt_path = which('main_cva.m');
% parrent_path = strrep(crt_path, 'Code\main_cva.m', '');
% data_file_path = strcat(parrent_path, 'Data\');
% 
% % read data from directory
% swap_rate_raw = xlsread(strcat(data_file_path, 'us_swap_rate.xlsx'));
% swap_tenor = swap_rate_raw(:,1);
% swap_rates = swap_rate_raw(:,2);

% Settle = datenum('30-June-2016');
% 
Tenor = 12*swap_tenor;
ZeroRates = swap_rates;

Settle = datenum('30-June-2016');

% Tenor = [3 6 12 5*12 7*12 10*12 20*12 30*12]';
% ZeroRates = [0.033 0.034 0.035 0.040 0.042 0.044 0.048 0.0475]';

ZeroDates = datemnth(Settle,Tenor);
Compounding = 2;
Basis = 0;
RateSpec = intenvset('StartDates', Settle,'EndDates', ZeroDates, ...
    'Rates', ZeroRates,'Compounding',Compounding,'Basis',Basis);

figure;
plot(ZeroDates, ZeroRates, 'o-');
xlabel('Date');
datetick('keeplimits');
ylabel('Zero rate');
grid on;
title('Yield Curve at Settle Date');

% Set Changeable Simulation Parameters
% Number of Monte Carlo simulations
numScenarios = 1000;

% Compute monthly simulation dates, then quarterly dates later.
simulationDates = datemnth(Settle,0:12);
simulationDates = [simulationDates datemnth(simulationDates(end),3:3:74)]';
numDates = numel(simulationDates);

% Compute Floating Reset Dates
floatDates = cfdates(Settle-360,swaps.Maturity,swaps.Period);
swaps.FloatingResetDates = zeros(numSwaps,numDates);
for i = numDates:-1:1
    thisDate = simulationDates(i);
    floatDates(floatDates > thisDate) = 0;
    swaps.FloatingResetDates(:,i) = max(floatDates,[],2);
end

% Setup Hull-White Single Factor Model
Alpha = 0.2;
Sigma = 0.015;

hw1 = HullWhite1F(RateSpec,Alpha,Sigma);

% Simulate Scenarios
% Use reproducible random number generator (vary the seed to produce
% different random scenarios).
prevRNG = rng(0, 'twister');

dt = diff(yearfrac(Settle,simulationDates,1));
nPeriods = numel(dt);
scenarios = hw1.simTermStructs(nPeriods, ...
    'nTrials',numScenarios, ...
    'deltaTime',dt);

% Restore random number generator state
rng(prevRNG);

% Compute the discount factors through each realized interest rate
% scenario.
dfactors = ones(numDates,numScenarios);
for i = 2:numDates
    tenorDates = datemnth(simulationDates(i-1),Tenor);
    rateAtNextSimDate = interp1(tenorDates,squeeze(scenarios(i-1,:,:)), ...
        simulationDates(i),'linear','extrap');
    % Compute D(t1,t2)
    dfactors(i,:) = zero2disc(rateAtNextSimDate, ...
        repmat(simulationDates(i),1,numScenarios),simulationDates(i-1),-1,3);
end
dfactors = cumprod(dfactors,1);

% Inspect a Scenario
i = 20;
figure;
surf(Tenor, simulationDates, scenarios(:,:,i))
axis tight
datetick('y','mmmyy');
xlabel('Tenor (Months)');
ylabel('Observation Date');
zlabel('Rates');
ax = gca;
ax.View = [-49 32];
title(sprintf('Scenario %d Yield Curve Evolution\n',i));

% Compute Mark to Market Swap Prices
% Compute all mark-to-market values for this scenario.  We use an
% approximation function here to improve performance.
values = hcomputeMTMValues(swaps,simulationDates,scenarios,Tenor);


% Inspect Scenario Prices
i = 32;
figure;
plot(simulationDates, values(:,:,i));
datetick;
ylabel('Mark-To-Market Price');
title(sprintf('Swap prices along scenario %d', i));

% Visualize Simulated Portfolio Values
% View portfolio value over time
figure;
totalPortValues = squeeze(sum(values(:,end,:), 2));
plot(simulationDates,totalPortValues);
title('Total MTM Portfolio Value for All Scenarios');
datetick('x','mmmyy')
legend
ylabel('Portfolio Value ($)')
xlabel('Simulation Dates')

% Compute Exposure by Counterparty
[exposures, expcpty] = creditexposures(values,swaps.Counterparty, ...
    'NettingID',swaps.NettingID);
% View portfolio exposure over time
figure;
totalPortExposure = squeeze(sum(exposures,2));
plot(simulationDates,totalPortExposure);
title('Portfolio Exposure for All Scenarios');
datetick('x','mmmyy')
ylabel('Exposure ($)')
xlabel('Simulation Dates')

% Exposure Profiles
% Compute entire portfolio exposure
portExposures = sum(exposures,2);

% Compute exposure profiles for each counterparty and entire portfolio
cpProfiles = exposureprofiles(simulationDates,exposures);
portProfiles = exposureprofiles(simulationDates,portExposures);

% Visualize portfolio exposure profiles
figure;
plot(simulationDates,portProfiles.PFE, ...
    simulationDates,portProfiles.MPFE * ones(numDates,1), ...
    simulationDates,portProfiles.EE, ...
    simulationDates,portProfiles.EPE * ones(numDates,1), ...
    simulationDates,portProfiles.EffEE, ...
    simulationDates,portProfiles.EffEPE * ones(numDates,1));
legend({'PFE (95%)','Max PFE','Exp Exposure (EE)','Time-Avg EE (EPE)', ...
    'Max past EE (EffEE)','Time-Avg EffEE (EffEPE)'})

datetick('x','mmmyy')
title('Portfolio Exposure Profiles');
ylabel('Exposure ($)')
xlabel('Simulation Dates')

cpIdx = find(expcpty == 6);
figure;
plot(simulationDates,cpProfiles(cpIdx).PFE, ...
    simulationDates,cpProfiles(cpIdx).MPFE * ones(numDates,1), ...
    simulationDates,cpProfiles(cpIdx).EE, ...
    simulationDates,cpProfiles(cpIdx).EPE * ones(numDates,1), ...
    simulationDates,cpProfiles(cpIdx).EffEE, ...
    simulationDates,cpProfiles(cpIdx).EffEPE * ones(numDates,1));
legend({'PFE (95%)','Max PFE','Exp Exposure (EE)','Time-Avg EE (EPE)', ...
    'Max past EE (EffEE)','Time-Avg EffEE (EffEPE)'})

datetick('x','mmmyy','keeplimits')
title(sprintf('Counterparty %d Exposure Profiles',cpIdx));
ylabel('Exposure ($)')
xlabel('Simulation Dates')

% Discounted Exposures
% Get discounted exposures per counterparty, for each scenario
discExp = zeros(size(exposures));
for i = 1:numScenarios
    discExp(:,:,i) = bsxfun(@times,dfactors(:,i),exposures(:,:,i));
end

% Discounted expected exposure
discProfiles = exposureprofiles(simulationDates,discExp, ...
    'ProfileSpec','EE');

% Aggregate the discounted EE for each counterparty into a matrix
discEE = [discProfiles.EE];

% Portfolio discounted EE
figure;
plot(simulationDates,sum(discEE,2))
datetick('x','mmmyy','keeplimits')
title('Discounted Expected Exposure for Portfolio');
ylabel('Discounted Exposure ($)')
xlabel('Simulation Dates')

% Counterparty discounted EE
figure;
plot(simulationDates,discEE)
datetick('x','mmmyy','keeplimits')
title('Discounted Expected Exposure for Each Counterparty');
ylabel('Discounted Exposure ($)')
xlabel('Simulation Dates')
legend('Counterparty 1', 'Counterparty 2', 'Counterparty 3', 'Counterparty 4','Counterparty 5','Counterparty 6')

% Calibrating Probability of Default Curve for Each Counterparty
% Import CDS market information for each counterparty
CDS = readtable(strcat(data_file_path, swapFile),'Sheet','CDS Spreads');
disp(CDS);
CDSDates = datenum(CDS.Date);
CDSSpreads = table2array(CDS(:,2:end));

ZeroData = [RateSpec.EndDates RateSpec.Rates];

% Calibrate default probabilities for each counterparty
DefProb = zeros(length(simulationDates), size(CDSSpreads,2));
for i = 1:size(DefProb,2)
    probData = cdsbootstrap(ZeroData, [CDSDates CDSSpreads(:,i)], ...
        Settle, 'probDates', simulationDates);
    DefProb(:,i) = probData(:,2);
end

% We plot of the cumulative probability of default for each counterparty.
figure;
plot(simulationDates,DefProb)
title('Default Probability Curve for Each Counterparty');
xlabel('Date');
grid on;
ylabel('Cumulative Probability')
datetick('x','mmmyy')
ylabel('Probability of Default')
xlabel('Simulation Dates')

% CVA Computation
Recovery = 0.4;
CVA = (1-Recovery) * sum(discEE(2:end,:) .* diff(DefProb));
for i = 1:numel(CVA)
    fprintf('CVA for counterparty %d = $%.2f\n',i,CVA(i));
end

figure;
bar(CVA);
title('CVA for each counterparty');
xlabel('Counterparty');
ylabel('CVA $');
grid on;