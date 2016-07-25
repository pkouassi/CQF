function outVol = blackvolbyrebonato(inCurve,VolFunc,CorrMat,ExerciseDate,Maturity,varargin)
%BLACKVOLBYREBONATO Compute black volatility for LMM using Rebonato formula
%
% Syntax:
%
%   outVol = blackvolbyrebonato(inCurve,VolFun,CorrMat,ExerciseDate,Maturity)
%   outVol = blackvolbyrebonato(inCurve,VolFun,CorrMat,ExerciseDate,Maturity,...
%                            'name1','val1')
%
% Description:
%
%   Compute black volatility for a swaption with the input exercise date
%   and maturity for a Libor Market Model specified by the cell array
%   of volatility function handles, and a correlation matrix for the
%   following Libor Market Model -- where the drift is
%   evolved using the spot measure.
%
%   dF_i/F_i = u_i + vol*dW_i, where dW is a i dimensional Brownian motion
%
%   and the drift terms are computed using the Spot measure.
%
% Inputs:
%
%   ZeroCurve - IRDataCurve or RateSpec. This is the zero curve that is
%               used to evolve the path of future interest rates.
%
%   VolFunc - NumRates-1 X 1 cell array of function handles. Each function
%             handle must take time as an input and return a scalar
%             volatility.
%
%   Correlation - NumRates-1 X NumRates-1 correlation matrix
%
%   ExerciseDate - NumSwaptions x 1 vector of serial date numbers or date strings
%                  containing the swaption exercise dates.
%
%   Maturity - NumSwaptions x 1 vector of serial date numbers or date strings
%              containing the swap maturity dates.
%
% Optional Inputs:
%
%   Period - Scalar, compounding frequency of curve and reset of swaptions
%            -- default is 2.
%
% Example:
%
%   Settle = datenum('11-Aug-2004');
% 
%   % Zero Curve
%   CurveTimes = (1:10)';
%   CurveDates = daysadd(Settle,360*CurveTimes,1);
% 
%   ZeroRates = [0.03 0.033 0.036 0.038 0.04 0.042 0.043 0.044 0.045 0.046]';
% 
%   % Construct an IRCurve
%   irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates);
% 
%   LMMVolFunc = @(a,t) (a(1)*t + a(2)).*exp(-a(3)*t) + a(4);
%   LMMVolParams = [.3 -.02 .7 .14];
% 
%   numRates = length(ZeroRates);
%   VolFunc(1:numRates-1) = {@(t) LMMVolFunc(LMMVolParams,t)};
% 
% 	Beta = .08;
%   CorrFunc = @(i,j,Beta) exp(-Beta*abs(i-j));
%   CorrMat = CorrFunc(meshgrid(1:numRates-1)',meshgrid(1:numRates-1),Beta);
% 
%   ExerciseDate = datenum('11-Aug-2009');
%   Maturity = daysadd(ExerciseDate,360*[3;4],1);
% 
%   Vol = blackvolbyrebonato(irdc,VolFunc,CorrMat,ExerciseDate,Maturity,'Period',1)
%
% Reference:
%
%   [1] Brigo, D and F. Mercurio. Interest Rate Models - Theory and
%   Practice. Springer Finance, 2006.
%
% See also LIBORMARKETMODEL

% Copyright 1999-2012 The MathWorks, Inc.

narginchk(5, 7);

if isafin(inCurve,'RateSpec')
    inCurve = IRDataCurve('Zero',inCurve.ValuationDate,inCurve.EndDates,...
        inCurve.Rates,'Basis',inCurve.Basis,'Compounding',inCurve.Compounding);
elseif ~isa(inCurve,'IRDataCurve')
    error(message('fininst:blackvolbyrebonato:invalidCurve'));
end

if iscell(VolFunc)
    if ~all(cellfun(@(x) isa(x,'function_handle'),VolFunc))
        error(message('fininst:blackvolbyrebonato:invalidVolFunc'));
    end
else
    error(message('fininst:blackvolbyrebonato:invalidVolFunc'));
end

if ~isnumeric(CorrMat) 
   error(message('fininst:blackvolbyrebonato:invalidCorrelation'));
end

ExerciseDate = datenum(ExerciseDate);
Maturity = datenum(Maturity);

if any(ExerciseDate > Maturity)
   error(message('fininst:blackvolbyrebonato:MaturityBeforeSettle'));
end

p = inputParser;
p.addParamValue('period',2);

try
    p.parse(varargin{:});
catch ME
    newMsg = message('fininst:blackvolbyrebonato:optionalInputError');
    newME = MException(newMsg.Identifier, getString(newMsg));
    newME = addCause(newME, ME);
    throw(newME)
end

Period = p.Results.period;

if length(VolFunc) ~= length(CorrMat)
    error(message('fininst:blackvolbyrebonato:matchingVolFuncCorrMat'));
end

try
    [ExerciseDate, Maturity] = finargsz(1, ExerciseDate(:),Maturity(:));
catch ME
    throwAsCaller(ME)
end

nRates = length(VolFunc)+1;

TExercise = round(yearfrac(inCurve.Settle,ExerciseDate,inCurve.Basis));
Tenor = round(yearfrac(ExerciseDate,Maturity,inCurve.Basis));

MaxSwapMat = round(yearfrac(inCurve.Settle,max(Maturity),inCurve.Basis));

if MaxSwapMat*Period > (nRates+1)
    error(message('fininst:blackvolbyrebonato:matchingVolFuncMaturity'));
end

CFDates = daysadd(inCurve.Settle,360*(1:(1/Period):nRates*Period),1);
ForwardRates = inCurve.getForwardRates(CFDates);
ForwardTimes = yearfrac(inCurve.Settle,CFDates,inCurve.Basis);

delta = [ForwardTimes(1);diff(ForwardTimes)];

wtemp = delta.*cumprod(1./(1 + ForwardRates.*delta));

% Compute for each swaption maturity
nSwaptions = length(ExerciseDate);
outVol = zeros(nSwaptions,1);

for iSwaption=1:nSwaptions
    
    alpha = TExercise(iSwaption)*Period;
    beta = TExercise(iSwaption)*Period + Tenor(iSwaption)*Period;
    
    % Compute weights
    w = wtemp(alpha+1:beta)./sum(wtemp(alpha+1:beta));
    
    SwapRate = sum(ForwardRates(alpha+1:beta).*w);
    SumTerm = 0;
    for i=alpha+1:beta
        for j=alpha+1:beta
            VolTerm = integral(@(x) VolFunc{i-1}(ForwardTimes(i-1)-x).*VolFunc{j-1}(ForwardTimes(j-1)-x),0,TExercise(iSwaption));
            SumTerm = SumTerm + (w(i-alpha)*w(j-alpha)*ForwardRates(i-1)*ForwardRates(j-1)*CorrMat(i-1,j-1))./ ...
                                        SwapRate.^2*VolTerm;
        end
    end
    outVol(iSwaption) = sqrt(SumTerm/TExercise(iSwaption));
end