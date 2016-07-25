function swaption_price = bump_curve_pricer(Settle,SwapTenor,ExerciseDates,SwapRates,BumpSize,BlackVol, InstrumentStrike)

CurveDates = daysadd(Settle,360*(SwapTenor'),1);
ZeroRates = max((SwapRates+BumpSize)/100,0);
irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates);
RateSpec = ...
    intenvset('Rates',ZeroRates,'EndDates',CurveDates,'StartDate',Settle);

% subsample of the whole vol surface
subsmpl = [1 3 5 7 10 13 15];
SwapTenor = SwapTenor(subsmpl);

EurExDatesFull = repmat(daysadd(Settle,ExerciseDates*360,1)',...
    length(SwapTenor),1);
EurMatFull = reshape(daysadd(EurExDatesFull,...
    repmat(360*SwapTenor,1,length(ExerciseDates)),1),size(EurExDatesFull));


nRates = 60;

CorrFunc = @(i,j,Beta) exp(-Beta*abs(i-j));

objfun = @(x) BlackVol(:) - blackvolbyrebonato(RateSpec,...
    repmat({@(t) ones(size(t)).*(x(1)*t + x(2)).*exp(-x(3)*t) + x(4)},nRates-1,1),...
    CorrFunc(meshgrid(1:nRates-1)',meshgrid(1:nRates-1),x(5)),...
    EurExDatesFull(:),EurMatFull(:),'Period',1);

options = optimset('disp','iter','MaxFunEvals',1000,'TolFun',1e-5);

x0 = [.2 .05 1 .05 .2];
lb = [0 0 .5 0 .01];
ub = [1 1 2 .3 1];
LMMparams = lsqnonlin(objfun,x0,lb,ub,options);

a = LMMparams(1);
b = LMMparams(2);
c = LMMparams(3);
d = LMMparams(4);

Beta = LMMparams(5);
CorrelationMatrix = CorrFunc(meshgrid(1:nRates-1)',meshgrid(1:nRates-1),Beta);
VolFunc = repmat({@(t) ones(size(t)).*(a*t + b).*exp(-c*t) + d},nRates-1,1);

nPeriods = 10;
DeltaTime = 1;
nTrials = 200;

Tenor = (1:10)';

SimDates = daysadd(Settle,360*DeltaTime*(0:nPeriods),1);
SimTimes = diff(yearfrac(SimDates(1),SimDates));

exRow = 5;
endCol = 5;

LMM = LiborMarketModel(irdc,VolFunc,CorrelationMatrix,'Period',1);
[LMMZeroRates, ~] = LMM.simTermStructs(nPeriods,'nTrials',nTrials);

DF = exp(bsxfun(@times,-LMMZeroRates(:,1:10,:),repmat(Tenor',[nPeriods+1 1])));
SwapRate = (1 - DF(exRow,endCol,:))./sum(bsxfun(@times,1,DF(exRow,1:endCol,:)));
PayoffValue = 100*max(SwapRate-InstrumentStrike,0).*sum(bsxfun(@times,1,DF(exRow,1:endCol,:)));
RealizedDF = prod(exp(bsxfun(@times,-LMMZeroRates(2:exRow,1,:),SimTimes(1:exRow-1))),1);
swaption_price = mean(RealizedDF.*PayoffValue);

end

