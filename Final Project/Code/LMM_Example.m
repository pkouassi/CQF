Settle = datenum('15-Dec-2007');
CurveTimes = [1:5 7 10 20]';
ZeroRates = [.01 .018 .024 .029 .033 .034 .035 .034]';
CurveDates = daysadd(Settle,360*CurveTimes,1);

irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates);

LMMVolFunc = @(a,t) (a(1)*t + a(2)).*exp(-a(3)*t) + a(4);
LMMVolParams = [.3 -.02 .7 .14];
  
numRates = 20;
VolFunc(1:numRates-1) = {@(t) LMMVolFunc(LMMVolParams,t)};
  
Beta = .08;
CorrFunc = @(i,j,Beta) exp(-Beta*abs(i-j));
Correlation = CorrFunc(meshgrid(1:numRates-1)',meshgrid(1:numRates-1),Beta);
  
LMM = LiborMarketModel(irdc,VolFunc,Correlation,'Period',1);