Settle = datenum('21-Jul-2008');

% Zero Curve
CurveDates = daysadd(Settle,360*([1 3 5 7 10 20]),1);
ZeroRates = [1.9 2.6 3.1 3.5 4 4.3]'/100;

plot(CurveDates,ZeroRates)
datetick
title(['Zero Curve for ' datestr(Settle)]);

irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates);

RateSpec = intenvset('Rates',ZeroRates,'EndDates',CurveDates,'StartDate',Settle)

InstrumentExerciseDate = datenum('21-Jul-2013');
InstrumentMaturity = datenum('21-Jul-2018');
InstrumentStrike = .045;

SwaptionBlackVol = [22 21 19 17 15 13 12
    21 19 17 16 15 13 11
    20 18 16 15 14 12 11
    19 17 15 14 13 12 10
    18 16 14 13 12 11 10
    15 14 13 12 12 11 10
    13 13 12 11 11 10 9]/100;
ExerciseDates = [1:5 7 10];
Tenors = [1:5 7 10];

EurExDatesFull = repmat(daysadd(Settle,ExerciseDates*360,1)',...
    length(Tenors),1);
EurMatFull = reshape(daysadd(EurExDatesFull,...
    repmat(360*Tenors,1,length(ExerciseDates)),1),size(EurExDatesFull));

% Find the swaptions that expire on or before the maturity date of the
% sample swaption
relidx = find(EurMatFull <= InstrumentMaturity);


SwaptionBlackPrices = zeros(size(SwaptionBlackVol));
SwaptionStrike = zeros(size(SwaptionBlackVol));

for iSwaption=1:length(ExerciseDates)
    for iTenor=1:length(Tenors)
        [~,SwaptionStrike(iTenor,iSwaption)] = swapbyzero(RateSpec,[NaN 0], Settle, EurMatFull(iTenor,iSwaption),...
            'StartDate',EurExDatesFull(iTenor,iSwaption),'LegReset',[1 1]);
        SwaptionBlackPrices(iTenor,iSwaption) = swaptionbyblk(RateSpec, 'call', SwaptionStrike(iTenor,iSwaption),Settle, ...
            EurExDatesFull(iTenor,iSwaption), EurMatFull(iTenor,iSwaption), SwaptionBlackVol(iTenor,iSwaption));
    end
end

nPeriods = 5;
DeltaTime = 1;
nTrials = 1000;

Tenor = (1:10)';

SimDates = daysadd(Settle,360*DeltaTime*(0:nPeriods),1)
SimTimes = diff(yearfrac(SimDates(1),SimDates))

% For 1 year periods and an evenly spaced tenor, the exercise row will be
% the sixth row and the swaption maturity will be the 5th column
exRow = 6;
endCol = 5;


nRates = 10;

CorrFunc = @(i,j,Beta) exp(-Beta*abs(i-j));

objfun = @(x) SwaptionBlackVol(relidx) - blackvolbyrebonato(RateSpec,...
    repmat({@(t) ones(size(t)).*(x(1)*t + x(2)).*exp(-x(3)*t) + x(4)},nRates-1,1),...
    CorrFunc(meshgrid(1:nRates-1)',meshgrid(1:nRates-1),x(5)),...
    EurExDatesFull(relidx),EurMatFull(relidx),'Period',1);

options = optimset('disp','iter','MaxFunEvals',1000,'TolFun',1e-5);

x0 = [.2 .05 1 .05 .2];
lb = [0 0 .5 0 .01];
ub = [1 1 2 .3 1];
LMMparams = lsqnonlin(objfun,x0,lb,ub,options)

a = LMMparams(1);
b = LMMparams(2);
c = LMMparams(3);
d = LMMparams(4);

Beta = LMMparams(5);

VolFunc = repmat({@(t) ones(size(t)).*(a*t + b).*exp(-c*t) + d},nRates-1,1);

CorrelationMatrix = CorrFunc(meshgrid(1:nRates-1)',meshgrid(1:nRates-1),Beta);
LMM = LiborMarketModel(irdc,VolFunc,CorrelationMatrix,'Period',1);

[LMMZeroRates, ForwardRates] = LMM.simTermStructs(nPeriods,'nTrials',nTrials);

trialIdx = 1;
figure
tmpPlotData = LMMZeroRates(:,:,trialIdx);
tmpPlotData(tmpPlotData == 0) = NaN;
surf(Tenor,SimDates,tmpPlotData)
title(['Evolution of the Zero Curve for Trial:' num2str(trialIdx) ' of LIBOR Market Model'])
xlabel('Tenor (Years)')

DF = exp(bsxfun(@times,-LMMZeroRates,repmat(Tenor',[nPeriods+1 1])));
SwapRate = (1 - DF(exRow,endCol,:))./sum(bsxfun(@times,1,DF(exRow,1:endCol,:)));
PayoffValue = 100*max(SwapRate-InstrumentStrike,0).*sum(bsxfun(@times,1,DF(exRow,1:endCol,:)));
RealizedDF = prod(exp(bsxfun(@times,-LMMZeroRates(2:exRow,1,:),SimTimes(1:exRow-1))),1);
LMM_SwaptionPrice = mean(RealizedDF.*PayoffValue)