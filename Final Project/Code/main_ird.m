% main function for the final projection

%% Part I Interest Rate Derivatives

% define workplace
crt_path = which('main_ird.m');
parrent_path = strrep(crt_path, 'Code\main_ird.m', '');
data_file_path = strcat(parrent_path, 'Data\');

% read data from directory
swap_rate_raw = xlsread(strcat(data_file_path, 'us_swap_rate.xlsx'));
swap_tenor = swap_rate_raw(:,1);
swap_rates = swap_rate_raw(:,2);

% construct zero curve and forward curve
Settle = datenum('30-June-2016');
CurveDates = daysadd(Settle,360*(swap_tenor'),1);
ZeroRates = swap_rates/100;
[ForwardRates, ~] = zero2fwd(ZeroRates, CurveDates, Settle);

figure
subplot(1,2,1)
plot(CurveDates,ZeroRates,'LineWidth',2)
datetick
title(['Zero Curve for ' datestr(Settle)]);
subplot(1,2,2)
plot(CurveDates,ForwardRates,'LineWidth',2)
datetick
title(['Forward Curve for ' datestr(Settle)]);

irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates);
RateSpec = ...
    intenvset('Rates',ZeroRates,'EndDates',CurveDates,'StartDate',Settle);

% read in swaption data
swaption_normal_vol = ...
    xlsread(strcat(data_file_path, 'swaption_normal_vol.xlsx')) / 100;

ExerciseDates = [1:10 12 15:5:30];
Tenors = [1:10 12 15:5:30];
figure
surf(ExerciseDates, Tenors, swaption_normal_vol);
xlabel('Option Terms')
ylabel('Swap Terms')
title('Swaption Normal Volatilities as of 2016-06-30')

% subsample of the whole vol surface
subsmpl = [1 3 5 7 10 13 15];
swaption_normal_vol = swaption_normal_vol(subsmpl,subsmpl);
ExerciseDates = ExerciseDates(subsmpl);
Tenors = Tenors(subsmpl);


EurExDatesFull = repmat(daysadd(Settle,ExerciseDates*360,1)',...
    length(Tenors),1);
EurMatFull = reshape(daysadd(EurExDatesFull,...
    repmat(360*Tenors,1,length(ExerciseDates)),1),size(EurExDatesFull));

% swaption normal vol to black vol
total_swaption = size(swaption_normal_vol,1) * size(swaption_normal_vol, 2);
swaption_strike = zeros(total_swaption,1);
swaption_black_vol = zeros(total_swaption,1);
for swpind = 1:total_swaption
[~,swaption_strike(swpind)] = swapbyzero(RateSpec,[NaN 0], Settle, EurMatFull(swpind),...
                'StartDate',EurExDatesFull(swpind),'LegReset',[1 1]);
        swaption_black_vol(swpind) = ...
                (swaption_normal_vol(swpind)/ ...
                swaption_strike(swpind))/100; 
end

% calibrate to market 
nRates = 60;

CorrFunc = @(i,j,Beta) exp(-Beta*abs(i-j));

objfun = @(x) swaption_black_vol(:) - blackvolbyrebonato(RateSpec,...
    repmat({@(t) ones(size(t)).*(x(1)*t + x(2)).*exp(-x(3)*t) + x(4)},nRates-1,1),...
    CorrFunc(meshgrid(1:nRates-1)',meshgrid(1:nRates-1),x(5)),...
    EurExDatesFull(:),EurMatFull(:),'Period',1);

options = optimset('disp','iter','MaxFunEvals',1000,'TolFun',1e-5);

x0 = [.2 .05 1 .05 .2];
lb = [0 0 .5 0 .01];
ub = [1 1 2 .3 1];
LMMparams = lsqnonlin(objfun,x0,lb,ub,options);

% LMM pricer 
a = LMMparams(1);
b = LMMparams(2);
c = LMMparams(3);
d = LMMparams(4);

Beta = LMMparams(5);

VolFunc = repmat({@(t) ones(size(t)).*(a*t + b).*exp(-c*t) + d},nRates-1,1);

figure
fplot(VolFunc{2},[0 30])
title('Volatility Function from Calibrated LMM')

% simulation
nPeriods = 10;
DeltaTime = 1;
nTrials = 1000;

Tenor = (1:10)';

SimDates = daysadd(Settle,360*DeltaTime*(0:nPeriods),1);
SimTimes = diff(yearfrac(SimDates(1),SimDates));

exRow = 5;
endCol = 5;

% correlation matrix
CorrelationMatrix = CorrFunc(meshgrid(1:nRates-1)',meshgrid(1:nRates-1),Beta);
disp('Correlation Matrix');
fprintf([repmat('%1.3f ',1,length(CorrelationMatrix)) ' \n'],CorrelationMatrix)
LMM = LiborMarketModel(irdc,VolFunc,CorrelationMatrix,'Period',1);
[LMMZeroRates, ForwardRates] = LMM.simTermStructs(nPeriods,'nTrials',nTrials);

% output scenario: zero curves and forward curves
trialIdx = 1;
figure
subplot(1,2,1)
tmpPlotData = LMMZeroRates(:,:,trialIdx);
tmpPlotData(tmpPlotData == 0) = NaN;
surf(Tenor,SimDates,tmpPlotData(:,1:10))
datetick('y')
title(['Evolution of the Zero Curve for Trial:' num2str(trialIdx) ' of LIBOR Market Model'])
xlabel('Tenor (Years)')
ylabel('Simulation Dates')
subplot(1,2,2)
tmpPlotData = ForwardRates(:,:,trialIdx);
tmpPlotData(tmpPlotData == 0) = NaN;
surf(Tenor,SimDates,tmpPlotData(:,1:10))
datetick('y')
title(['Evolution of the Forward Curve for Trial:' num2str(trialIdx) ' of LIBOR Market Model'])
xlabel('Tenor (Years)')
ylabel('Simulation Dates')

% price of the swaption
InstrumentStrike = .025;
DF = exp(bsxfun(@times,-LMMZeroRates(:,1:10,:),repmat(Tenor',[nPeriods+1 1])));
SwapRate = (1 - DF(exRow,endCol,:))./sum(bsxfun(@times,1,DF(exRow,1:endCol,:)));
PayoffValue = 100*max(SwapRate-InstrumentStrike,0).*sum(bsxfun(@times,1,DF(exRow,1:endCol,:)));
RealizedDF = prod(exp(bsxfun(@times,-LMMZeroRates(2:exRow,1,:),SimTimes(1:exRow-1))),1);
LMM_SwaptionPrice = mean(RealizedDF.*PayoffValue);

% Bermudans option pricing
% LMMTenorBer = (1:11)';
% BermudanExerciseDates = daysadd(Settle,360*(1:11),1);
% BermudanMaturity = datenum('30-June-2021');
% BermudanStrike = .045;
% LMMBermPrice = hBermudanSwaption(LMMZeroRates(6:11,1:11,:),SimDates(6:11),LMMTenorBer,.025,SimDates(1:5),BermudanMaturity);


% market risk analysis
swaption_price_0 = ...
    bump_curve_pricer(Settle,swap_tenor,ExerciseDates,swap_rates,0,swaption_black_vol, .025);
swaption_price_300 = ...
    bump_curve_pricer(Settle,swap_tenor,ExerciseDates,swap_rates,3,swaption_black_vol, .025);
swaption_price_100 = ...
    bump_curve_pricer(Settle,swap_tenor,ExerciseDates,swap_rates,1,swaption_black_vol, .025);
swaption_price_50 = ...
    bump_curve_pricer(Settle,swap_tenor,ExerciseDates,swap_rates,0.5,swaption_black_vol, .025);
swaption_price_25 = ...
    bump_curve_pricer(Settle,swap_tenor,ExerciseDates,swap_rates,0.25,swaption_black_vol, .025);
swaption_price_10 = ...
    bump_curve_pricer(Settle,swap_tenor,ExerciseDates,swap_rates,0.1,swaption_black_vol, .025);
swaption_price_n10 = ...
    bump_curve_pricer(Settle,swap_tenor,ExerciseDates,swap_rates,-0.1,swaption_black_vol, .025);
swaption_price_n25 = ...
    bump_curve_pricer(Settle,swap_tenor,ExerciseDates,swap_rates,-0.25,swaption_black_vol, .025);
swaption_price_n50 = ...
    bump_curve_pricer(Settle,swap_tenor,ExerciseDates,swap_rates,-0.5,swaption_black_vol, .025);
% swaption_price_n100 = ...
%     bump_curve_pricer(Settle,swap_tenor,ExerciseDates,swap_rates,-1,swaption_black_vol, .025);

rho_swaption = (swaption_price_25 - swaption_price_n25)/50;
convexity_swaption = (swaption_price_25 - 2*swaption_price_0 + swaption_price_n25)/625;
shocks = [-50, -25, -10, 0, 10, 25, 50, 100, 300];
swaption_prices = [swaption_price_n50, swaption_price_n25, swaption_price_n10, ...
    swaption_price_0, swaption_price_10, swaption_price_25,swaption_price_50, ...
    swaption_price_100, swaption_price_300];
figure
plot(shocks, swaption_prices,'LineWidth',2);
title('Price Sensitivity of Interest Rate wrt Swaption Price');
