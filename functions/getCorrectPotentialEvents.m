function trialsIndices = getCorrectPotentialEvents(deltaFStacks, eventsPoints, xbegin, xresp, pseudoxresp)
%GETCORRECTPOTENTIALEVENTS Summary of this function goes here
%   Detailed explanation goes here
    
    % Iterate over Potential events traces to get the Critical Chi2 Value
    % from baseline pseudo-fit
    chi2Samples = [];
    for t = 1 : size(eventsPoints{1}, 1)
        trace = reshape(deltaFStacks{1}(eventsPoints{1}(t,3), eventsPoints{1}(t,2), :),1,[]);
        baseline = trace(xbegin-1:xresp-2) ;
        [fitresult, ~] = fitSignalTemplate(baseline, 1, pseudoxresp);
        % Chi square calculation
        fit = fitresult(1:length(baseline))';
        chi2Samples = [chi2Samples, sum((baseline - fit).^2 / fit.^2)];
    end
    alpha = 0.1;
    chi2Crit = prctile(chi2Samples, 100*(1-alpha));
    
    % Iterate over Potential events traces
    trialsIndices = cell(1);
    parfor T = 1 : length(eventsPoints)
        trialIndex = zeros(1, size(eventsPoints{T}, 1));
        for t = 1 : size(eventsPoints{T}, 1)
            trace = reshape(deltaFStacks{T}(eventsPoints{T}(t,3), eventsPoints{T}(t,2), :),1,[]);
                         
            % ===== Fit traces with Signal Template =====
            [fitresult, ~] = fitSignalTemplate(trace, xbegin, xresp);

            % ===== Check curves similarity =====
            % Chi square calculation
            fit = fitresult(1:length(trace))';
            chi2 = sum((trace(xresp:end) - fit(xresp:end)).^2 / var(fit(xresp:end)));
            
            if chi2 <= chi2Crit
                trialIndex(t) = 1;
            end
        end
        trialsIndices{T} = trialIndex;
    end
end

function [fitresult, gof] = fitSignalTemplate(trace, xbegin, xresp)
    x = (1:length(trace));
    [xData, yData] = prepareCurveData(x, double(trace)); 
    
    A = median(trace(xbegin:xresp-3));
    B = mean(trace(xresp-1: xresp+2)) / median(trace(xresp-4:xresp-2)) - 1;
    tau_decay = 7;
    
    eq =  'A * ((1 + (B * exp (-(x-xresp) ./ tau_decay))) .* (x>=xresp) + (x<xresp))';
    ft = fittype(eq, 'independent', 'x', 'dependent', 'y');

    excludedVals = false(1,length(trace));
    excludedVals(1:xbegin+1) = true;
    outliers = excludedata(xData,yData,'indices',excludedVals);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', ...
        [0 -2^16 1 xresp], 'Upper',[2^16 2^16 10 xresp],...
        'MaxIter', 3000, 'Exclude', outliers, 'Display', 'notify');
    opts.StartPoint = [A B tau_decay xresp];
    
    [fitresult, gof] = fit(xData, yData, ft, opts);
end
