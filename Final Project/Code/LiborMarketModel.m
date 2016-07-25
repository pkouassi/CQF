classdef LiborMarketModel
%LIBORMARKETMODEL Create a Libor Market Model
%
% Syntax:
%
%   OBJ = LiborMarketModel(ZeroCurve,VolFunc,Correlation)
%
% Description:
%
%   Create a Libor Market Model by specifying a zero curve, a cell array
%   of volatility function handles, and a correlation matrix for the
%   following Libor Market Model -- where the drift is
%   evolved using the spot measure.
%
%   dF_i/F_i = u_i + vol*dW_i, where dW is a i dimensional Brownian motion
%
%   and the drift terms are computed using the Spot measure.
%
% Properties:
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
%   NumFactors - Number of Brownian factors. Default is NaN, in which case
%                the number of factors is equal to the number of rates.
%
%   Period - Period of the forward rates. Default is 2, meaning forward
%            rates are spaced at 0, .5, 1, 1.5 etc.
%
%   Note that Correlation and VolFunc are sized with NumRates-1 since the
%   first rate is locked in and essentialy dead.
%
% Methods:
%
%   [ZeroRates, ForwardRates] = simTermStructs(nPeriods)
%   [ZeroRates, ForwardRates] = simTermStructs(nPeriods,'name1','val1')
%
% Example:
%
%   Settle = datenum('15-Dec-2007');
%   CurveTimes = [1:5 7 10 20]';
%   ZeroRates = [.01 .018 .024 .029 .033 .034 .035 .034]';
%   CurveDates = daysadd(Settle,360*CurveTimes,1);
%
%   irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates);
%
%   LMMVolFunc = @(a,t) (a(1)*t + a(2)).*exp(-a(3)*t) + a(4);
%   LMMVolParams = [.3 -.02 .7 .14];
% 
%   numRates = 20;
%   VolFunc(1:numRates-1) = {@(t) LMMVolFunc(LMMVolParams,t)};
% 
%   Beta = .08;
%   CorrFunc = @(i,j,Beta) exp(-Beta*abs(i-j));
%   Correlation = CorrFunc(meshgrid(1:numRates-1)',meshgrid(1:numRates-1),Beta);
% 
%   LMM = LiborMarketModel(irdc,VolFunc,Correlation,'Period',1);
% 
%   [ZeroRates, ForwardRates] = LMM.simTermStructs(10,'nTrials',100);
%
% Reference:
%
%   [1] Brigo, D and F. Mercurio. Interest Rate Models - Theory and
%   Practice. Springer Finance, 2006.
% 
% See also HULLWHITE1F, LINEARGAUSSIAN2F

% Copyright 1999-2012 The MathWorks, Inc.   
    
    properties
        ZeroCurve
        VolFunctions
        Correlation
        NumFactors
        Period
    end
    
    properties (Access = private)
        SDE
        Tenor
        ForwardRates
    end
    
    methods (Access = public)
        function obj = LiborMarketModel(inCurve,inVol,inCorr,varargin)
        %LIBORMARKETMODEL Create a Libor Market Model
            
            narginchk(3,7);
            
            % If RateSpec, convert to be an IRDataCurve
            if isafin(inCurve,'RateSpec')
                obj.ZeroCurve = IRDataCurve('Zero',inCurve.ValuationDate,inCurve.EndDates,...
                    inCurve.Rates,'Basis',inCurve.Basis,'Compounding',inCurve.Compounding);
            elseif isa(inCurve,'IRDataCurve')
                obj.ZeroCurve = inCurve;
            else
                error(message('fininst:LiborMarketModel:invalidCurve'));
            end
            
            if iscell(inVol)
                if all(cellfun(@(x) isa(x,'function_handle'),inVol))
                    obj.VolFunctions = inVol;
                else
                    error(message('fininst:LiborMarketModel:invalidVolFunc'));
                end
            else
                error(message('fininst:LiborMarketModel:invalidVolFunc'));
            end
            
            if isnumeric(inCorr)
                obj.Correlation = inCorr;
            else
                error(message('fininst:LiborMarketModel:invalidCorrelation'));
            end
            
            if length(obj.VolFunctions) ~= length(obj.Correlation)
                error(message('fininst:LiborMarketModel:matchingVolFuncCorrMat'));
            end
            
            p = inputParser;
            
            p.addParamValue('numfactors',NaN,@isscalar);
            p.addParamValue('period',2);
            
            try
                p.parse(varargin{:});
            catch ME
                newMsg = message('fininst:LiborMarketModel:optionalInputError');
                newME = MException(newMsg.Identifier, getString(newMsg));
                newME = addCause(newME, ME);
                throw(newME)
            end
            
            if isnan(p.Results.numfactors) || mod(p.Results.numfactors,1) == 0
                obj.NumFactors = p.Results.numfactors;
            else
                error(message('fininst:LiborMarketModel:invalidNumFactors'));
            end
            obj.Period = p.Results.period;
            
            nRates = length(obj.VolFunctions)+1;
            obj.Tenor = (1/obj.Period):(1/obj.Period):(nRates/obj.Period);
            
            TenorDates = daysadd(obj.ZeroCurve.Settle,360*obj.Tenor,1);
            obj.ForwardRates = obj.ZeroCurve.getForwardRates(TenorDates);
            
            ForwardTimes = yearfrac(obj.ZeroCurve.Settle,TenorDates,obj.ZeroCurve.Basis);
            
            % If no factors specified, construct the GBM object using the
            % correlation matrix specified.
            if isnan(obj.NumFactors)
                try
                    obj.SDE = gbm(@(t,L) DriftSpot(t,L,obj.VolFunctions,obj.Correlation,obj.Tenor),...
                    @(t,L) Diffusion(t,obj.VolFunctions,ForwardTimes),...
                    'Correlation',inCorr,'StartState',obj.ForwardRates(2:end));
                catch ME
                    error(message('fininst:LiborMarketModel:invalidVolFunc'));
                end
            else
                % Decompose the correlation matrix and hit the volatilities
                % with it
                [~,S,V] = svd(obj.Correlation);
                B = -(V*sqrt(S(:,1:obj.NumFactors)));
                try                    
                    obj.SDE = gbm(@(t,L) DriftSpot(t,L,obj.VolFunctions,obj.Correlation,obj.Tenor),...
                        @(t,L) Diffusion(0,obj.VolFunctions,ForwardTimes)*B,...
                        'StartState',obj.ForwardRates(2:end));
                catch ME
                    error(message('fininst:LiborMarketModel:invalidVolFunc'));
                end
            end
            
        end
        function [ZeroRates, ForwardRates] = simTermStructs(obj,nPeriods,varargin)
        %SIMTERMSTRUCTS Simulate Term Structures
        %
        % Syntax:
        %
        %   [ZeroRates, ForwardRates] = simTermStructs(nPeriods)
        %   [ZeroRates, ForwardRates] = simTermStructs(nPeriods,'name1','val1')
        %
        % Description:
        %
        %   Simulate future zero curve paths using the specified Libor
        %   Market Model.
        %
        % Required Input Arguments:
        %
        %   nPeriods - Number of simulation periods
        %
        % Option Input Arguments:
        %
        %   deltaTime - scalar time step betwen periods. Default is 1.
        %
        %   nTrials - scalar number of trials. Default is 1.
        % 
        %   antithetic - Boolean scalar flag indicating whether antithetic
        %                sampling is used to generate the Gaussian random
        %                variates that drive the zero-drift, unit-variance
        %                rate Brownian vector dW(t). See
        %                gbm/simBySolution for more information.
        %
        %   Z - Direct specification of the dependent random noise process
        %       used to generate the zero-drift, unit-variance rate Brownian
        %       vector dW(t) that drives the simulation. See
        %       gbm/simBySolution for more information.
        %
        %   Tenor - numeric vector of maturities to be computed at each
        %           time step. Default is the number of rates in the
        %           object as specified by the correlation and volatility
        %           functions.
        % 
        % Output Arguments:
        %
        %   ZeroRates - nPeriods X nTenors X nTrials matrix of simulated
        %               zero rate term structures.
        %
        %   ForwardRates - nPeriods X nTenors X nTrials matrix of simulated
        %               forward rate term structures.
        %
        % Example:
        %
        %   CurveTimes = [1:5 7 10 20]';
        %   ZeroRates = [.01 .018 .024 .029 .033 .034 .035 .034]';
        %
        %   LMMVolFunc = @(a,t) (a(1)*t + a(2)).*exp(-a(3)*t) + a(4);
        %   LMMVolParams = [.3 -.02 .7 .14];
        % 
        %   numForwardRates = length(ForwardRates);
        %   VolFunc(1:numForwardRates-1) = {@(t) LMMVolFunc(LMMVolParams,t)};
        % 
        %   Beta = .08;
        %   CorrFunc = @(i,j,Beta) exp(-Beta*abs(i-j));
        %   Correlation = CorrFunc(meshgrid(1:numForwardRates-1)',meshgrid(1:numForwardRates-1),Beta);
        % 
        %   LMM = LiborMarketModel([CurveTimes ForwardRates],VolFunc,Correlation);
        % 
        %   [ZeroRates, ForwardRates] = LMM.simTermStructs(20,'nTrials',100);
        
            p = inputParser;
            
            p.addParamValue('ntrials',1);
            p.addParamValue('deltatime',1);
            p.addParamValue('tenor',[]);
            p.addParamValue('antithetic',false);
            p.addParamValue('Z',[]);
            
            try
                p.parse(varargin{:});
            catch ME
                newMsg = message('fininst:LiborMarketModel:optionalInputError');
                newME = MException(newMsg.Identifier, getString(newMsg));
                newME = addCause(newME, ME);
                throw(newME)
            end
            
            nTrials = p.Results.ntrials;
            deltaTime = p.Results.deltatime;
            inTenor = p.Results.tenor;
            Antithetic = p.Results.antithetic;
            Z = p.Results.Z;
            
            % Simulate
            Paths = obj.SDE.simBySolution(nPeriods,'NTRIALS',nTrials,'DeltaTime',deltaTime,...
                'antithetic',Antithetic,'Z',Z);
            
            FullPaths = [repmat([obj.ForwardRates(1);zeros(nPeriods,1)],[1 1 nTrials]) Paths];
            
            % Kill off the rates once the simulation has moved past them
            ForwardRates = bsxfun(@times,FullPaths,triu(ones([size(FullPaths,1) size(FullPaths,2)])));
            
            % Shift values in matrix to account for the convention that as
            % we move forward in periods the columns represent time to
            % maturity and not one particular rate. This is to be
            % consistent with other functions like HullWhite1F and
            % LinearGaussian2F
            TenorTimes = [obj.Tenor(1) diff(obj.Tenor)];
            
            ZeroRates = zeros(size(ForwardRates));
            ZeroRates(1,:,:) = bsxfun(@rdivide,cumsum(bsxfun(@times,TenorTimes,ForwardRates(1,:,:)),2),obj.Tenor);
            
            for iPeriod=2:nPeriods+1
                ZeroRates(iPeriod,:,:) = [bsxfun(@rdivide,cumsum(bsxfun( ...
                    @times,TenorTimes(iPeriod:end),ForwardRates(iPeriod,iPeriod:end,:)),2) ...
                    ,obj.Tenor(1:end-iPeriod+1)) zeros(1,iPeriod-1,nTrials)];
                ForwardRates(iPeriod,:,:) = [ForwardRates(iPeriod,iPeriod:end,:) zeros(1,iPeriod-1,nTrials)];
            end  
            
            if ~isempty(inTenor)
                tmpZeroRates = ZeroRates;
                tmpForwardRates = ForwardRates;
                ZeroRates = zeros(nPeriods,length(inTenor),nTrials);
                ForwardRates = ZeroRates;
                for iTrial=1:nTrials
                    for iPeriod=1:nPeriods
                        ZeroRates(iPeriod,:,iTrial) = interp1(obj.Tenor,tmpZeroRates(iPeriod,:,iTrial),inTenor);
                        ForwardRates(iPeriod,:,iTrial) = interp1(obj.Tenor,tmpForwardRates(iPeriod,:,iTrial),inTenor);
                    end
                end
            end
        end
    end
end

function outDiffusion = Diffusion(t,VolFunctions,ForwardTimes)

outDiffusion = zeros(length(VolFunctions));

for volidx=1:length(VolFunctions)
    outDiffusion(volidx,volidx) = VolFunctions{volidx}(ForwardTimes(volidx)-t);
end

end

function Drifts = DriftSpot(t,L,VolFunc,CorrMat,ForwardTimes)

tau = diff(ForwardTimes);

Drifts = zeros(size(L));

tidx = find(ForwardTimes == t);

if isempty(tidx),tidx = 1;end

for kdx=1:length(L)
    SumTerm = zeros(size(L));
    for jdx=tidx:kdx
        SumTerm(jdx) = (L(jdx)*tau(jdx).*CorrMat(jdx,kdx).*VolFunc{jdx}(ForwardTimes(jdx)-t))./(1 + tau(jdx)*L(jdx));
    end
    Drifts(kdx) = VolFunc{kdx}(ForwardTimes(kdx)-t)*sum(SumTerm);
end

Drifts = diag(Drifts);

end