classdef matchFactors < handle
    %UNTITLED2 This class takes in internal and external balance factor data structures for a full raneg of pressure
    % The class is used to determine pressures for each factor where the
    % factors are most balanced
    %   Detailed explanation goes here
    
    properties
        balInt = [];
        balExt = [];
        HRE1 = [];
        HRE2 = [];
        HRE3 = [];
        
        SK0 = [];
        SK1 = [];
        SK2 = [];
    end
    
    methods
        % constructor to instantiate object with internal and external
        % balance factor data structures
        function obj = matchFactors(balInt, balExt)
            if nargin == 0
                error('Must instantiate with balance factor data structures');
            end
            
            obj.balInt = balInt;
            obj.balExt = balExt;
        end
        
        function icp = findBalancePressure(obj, balFactor, bTCD3D, bThresh)
            % determines the balance pressure for a given balance factor
            % bThresh = 1 (crosses 1)
            % bThresh = 0 (max)
            % bThresh = -1 (min)
            
            numPressures = size(obj.balInt, 1);
            numIterations = size(obj.balInt, 2);

            ratioFactor = zeros(1, numPressures);
            
            for a = 1:numPressures
                ratioFactor1 = ones(1, numIterations)*NaN;
                for b = 1:numIterations
                    if bTCD3D
                        if obj.balInt{a,b}.snr == 0 || obj.balExt{a,b}.snr == 0
                            % no valid data for this pressure/ iteration
                        else
                            % process
                            ratioFactor1(b) = eval(['obj.' balFactor '{a,b};']);
                        end
                    else
                        if obj.balInt{a,b}.snr == 0 || obj.balExt{a,b}.snr == 0
                            % no valid data for this pressure/ iteration
                        else
                            % process
                            ratioFactor1(b) = eval(['obj.balInt{a,b}.' balFactor '.mean / obj.balExt{a,b}.' balFactor '.mean;']);
                        end
                    end
                end
                ratioFactor(a) = nanmean(ratioFactor1);
                pressure(a) = obj.balInt{a,1}.pressure;
            end
            
            % fit 2nd order polynomial to factor
            pressureVec = min(pressure):max(pressure);
            idx = find(isnan(ratioFactor));
            ratioFactor(idx) = [];
            pressure(idx) = [];
            if isempty(pressure)
                error('Need to handle this scenario if it exists!!')
            end
            [P, S] = polyfit(pressure, ratioFactor, 2);
            polynomFit = polyval(P, pressureVec);
            
            % check for crossing of 1, max or min
            thresh = 1;
            switch bThresh
                case 1
                    thresh = 1;
                case -1
                    thresh = 0;
                case 0
                    %(swap max and min TODO)
            end
            
            if  max(polynomFit) >= thresh && min(polynomFit) <= thresh
                numCrossings = 0;
                icp = [];
                for a = 1:length(polynomFit)-1
                    if (polynomFit(a) <= thresh && polynomFit(a+1) >= thresh) || (polynomFit(a+1) <= thresh && polynomFit(a) >= thresh)
                        numCrossings = numCrossings + 1;
                        icp(numCrossings) = interp1( [polynomFit(a) polynomFit(a+1)], [pressureVec(a) pressureVec(a+1)], 1, 'linear');
                    end
                end
                icp = mean(icp);
            else
                %icp = NaN;
                [~, I] = max(polynomFit);
                if I == 1 || I == length(polynomFit)
                    [~, I] = min(polynomFit);
                end
                icp = pressureVec(I);
            end
        end
        
        function getBalanceTCD3D(obj)
            % get TCD3D balance factors
            numPressures = size(obj.balInt, 1);
            numIterations = size(obj.balInt, 2);
            
            for a = 1:numPressures
                for b = 1:numIterations
                    if obj.balInt{a,b}.snr == 0 || obj.balExt{a,b}.snr == 0
                        % no valid data for this pressure/ iteration
                    else
                        % process
                        obj.HRE1{a,b} = obj.balInt{a,b}.PI_depth1111.mean./obj.balExt{a,b}.PI_depth2222.mean;
                        obj.HRE2{a,b} = (sum((obj.balInt{a,b}.PI_depth1102.mean - obj.balExt{a,b}.PI_depth2202.mean).^2))^0.5;
                        obj.HRE3{a,b} = (sum((obj.balInt{a,b}.PI_depth1100.mean - obj.balExt{a,b}.PI_depth2200.mean).^2))^0.5;
                        
                        obj.SK1{a,b} = (sum(([obj.balExt{a,b}.H10ex.mean obj.balExt{a,b}.H21ex.mean obj.balExt{a,b}.H32ex.mean obj.balExt{a,b}.H43ex.mean] - [obj.balInt{a,b}.H10int.mean obj.balInt{a,b}.H21int.mean obj.balInt{a,b}.H32int.mean obj.balInt{a,b}.H43int.mean]).^2))^0.5;
                        obj.SK2{a,b} = (sum(([obj.balExt{a,b}.H10ex.mean obj.balExt{a,b}.H21ex.mean obj.balExt{a,b}.H32ex.mean] - [obj.balInt{a,b}.H10int.mean obj.balInt{a,b}.H21int.mean obj.balInt{a,b}.H32int.mean]).^2))^0.5;
                        obj.SK0{a,b} = (sum(([obj.balExt{a,b}.H21ex.mean obj.balExt{a,b}.H32ex.mean obj.balExt{a,b}.H43ex.mean] - [obj.balInt{a,b}.H21int.mean obj.balInt{a,b}.H32int.mean obj.balInt{a,b}.H43int.mean]).^2))^0.5;
                    end
                end
                
            end
        end
        
        function [wavDist, dPI] = compareShapes(obj)
            % This function synchronizes two balance factor objects and
            % normalizes each waveform and then computes and returns
            % distance metric between each waveform
            
            minAllowableDiff = 3;
            count = 1;
            wavDist = [];
            dPI = [];
            % Check to see if any available cycles and then sycnh to best
            % line up
            if obj.balInt.snr == 0 || obj.balExt.snr == 0
                % exit
            else
                numIntCycles = length(obj.balInt.cycles.PeriodStart);
                numExtCycles = length(obj.balExt.cycles.PeriodStart);
                for a = 1:numIntCycles
                    periodStartInt = obj.balInt.cycles.PeriodStart(a);
                    periodEndInt = obj.balInt.cycles.PeriodEnd(a);
                    % find closest match in ext
                    fStart = find(abs(obj.balExt.cycles.PeriodStart -  periodStartInt) <= minAllowableDiff);
                    fEnd = find(abs(obj.balExt.cycles.PeriodEnd -  periodEndInt) <= minAllowableDiff);
                    if ~isempty(fStart) && ~isempty(fEnd)
                        periodStart = max(periodStartInt, obj.balExt.cycles.PeriodStart(fStart(1))); 
                        periodEnd = min(periodEndInt, obj.balExt.cycles.PeriodEnd(fEnd(1))); 
                        
                        wav1 = obj.balInt.env(periodStart:periodEnd);
                        wav2 = obj.balExt.env(periodStart:periodEnd);
                        
                        wav1norm = wav1./max(wav1);
                        wav2norm = wav2./max(wav2);
                        
                        wavDist(count) = pdist([wav1norm; wav2norm], 'euclidean');
                        
                        pi1 = obj.balInt.getPI(wav1);
                        pi2 = obj.balInt.getPI(wav2);
                        
                        dPI(count) = pi1 - pi2;
                        count = count + 1;
                    end
                end
                
            end
        end
        
    end
end
    



            