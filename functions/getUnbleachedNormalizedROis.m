function [correctedROIs, normalizedROIs, fitCoeffs] = getUnbleachedNormalizedROis(ROIs, xbegin, xresp)
%GETUNBLEACHEDNORMALIZEDROIS Summary of this function goes here
    tic
    [correctedROIs, normalizedROIs, fitCoeffs] = deal(cell(length(ROIs),1));
    parfor N = 1 : size(ROIs, 1)
        [T, C] = size(ROIs{N});
        ROI = double(ROIs{N});

        [newROICorrected, newROINorm] = deal(zeros(T,C));
        [newROICorrected(:,1), newROINorm(:,1)] = deal(ROIs{N}(:,1));
        ROICoeffs = NaN(T, 9);

        for t = 1 : T
            trace = ROI(t, 2:end);

            % Linear fitting of baseline
            baseline = trace(xbegin : xresp-3);
            [lineFit, baselineDecayPerc] = getBaselineDecay(baseline, xbegin, xresp);

            % Negligible Bleaching
            if baselineDecayPerc <= 1
                % No Bleaching Correction
                [fitresult, pseudoFit] = fitSignalTemplate(trace, xbegin, xresp);
                newROICorrected(t,2:end) = trace;
                newROINorm(t,2:end) = normalizeTrace(trace, fitresult);
                ROICoeffs(t,:) = getFittingCoefficients(fitresult, pseudoFit, baselineDecayPerc, 'expResponse');
            % Small Bleaching
            elseif baselineDecayPerc > 1 && baselineDecayPerc <= 6
                % Linear-Exponential fitting of the trace
                [fitresult, pseudoFit] = fitLinearExponentialTrace(trace, xbegin, xresp, lineFit);
                [newROICorrected(t,2:end), newROINorm(t,2:end)] = getLineCorrectionAndNormalization(trace, fitresult, xbegin);
                ROICoeffs(t,:) = getFittingCoefficients(fitresult, pseudoFit, baselineDecayPerc, 'lineExp');
            % Relevant Bleaching    
            else % baselineDecayPerc > 6
                % Double Exponential fitting of the trace
                [fitresult, pseudoFit] = fitDoubleExponentialTrace(trace, xbegin, xresp);
                [newROICorrected(t,2:end), newROINorm(t,2:end)] = getExponentialCorrectionAndNormalization(trace, fitresult, xbegin);
                ROICoeffs(t,:) = getFittingCoefficients(fitresult, pseudoFit, baselineDecayPerc, 'DoubleExp');
            end
        end
        correctedROIs{N} = newROICorrected;
        normalizedROIs{N} = newROINorm;
        fitCoeffs{N} = ROICoeffs;
    end
    toc
    disp("All traces corrected")
end

function [lineFit, baselineDecayPerc] = getBaselineDecay(baseline, xbegin, xresp)
    baselineX = (xbegin: xresp-3);
    lineFit = fit(baselineX',baseline','poly1', 'Lower', [-inf 0], 'Upper',[0 inf]);
    fittedBaseline = lineFit(baselineX)';
    relatSlope = abs(fittedBaseline(end) - fittedBaseline(1));
    baselineDecayPerc = (relatSlope * 100) / median(fittedBaseline);
end

function [fitresult, pseudoFit] = fitSignalTemplate(trace, xbegin, xresp)
    x = (1:length(trace));
    [xData, yData] = prepareCurveData(x, double(trace)); 
    
    A = median(trace(xbegin:xresp-3));
    B = mean(trace(xresp: xresp+2)) / median(trace(xresp-4:xresp-2)) - 1;
    tau_decay = 7;
    
    eq =  '((1 + (B * exp (-(x-xresp) ./ tau_decay))) .* (x>=xresp) + (x<xresp)) * A';
    ft = fittype(eq, 'independent', 'x', 'dependent', 'y');

    excludedVals = false(1,length(trace));
    excludedVals(1:xbegin-1) = true;
    outliers = excludedata(xData,yData,'indices',excludedVals);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', ...
        [0 -1 2 xresp], 'Upper',[2^16 20 30 xresp],...
        'MaxIter', 3000, 'Exclude', outliers, 'Display', 'notify');
    
    opts.StartPoint = [A B tau_decay xresp];
    
    [fitresult, ~] = fit(xData, yData, ft, opts);

    % Placebo Fitting
    [xData, yData] = prepareCurveData(x(1:xresp-3), double(trace(1:xresp-3)));
    outliers = excludedata(xData,yData,'indices',excludedVals(1:xresp-3));
    coeffs = coeffvalues(fitresult);
    opts.Exclude = outliers;
    opts.Lower = [coeffs(1), -1, 2,  10];
    opts.Upper = [coeffs(1), 20, 30, 10];
    opts.StartPoint = [coeffs(1) B tau_decay 10];
    [pseudoFit, ~] = fit(xData, yData, ft, opts);
end

function [fitresult, pseudoFit] = fitLinearExponentialTrace(trace, xbegin, xresp, linearFit)
    linCoeffs = coeffvalues(linearFit);
    x = (1:length(trace));
    [xData, yData] = prepareCurveData(x, double(trace));

    M = linCoeffs(1);
    A = linCoeffs(2);
    B = mean(trace(xresp+1: xresp+3)) / mean(trace(xresp-4:xresp-2)) - 1;
    tau_decay = 7;
    
    eq =  '((1 + (B * exp (-(x-xresp) ./ tau_decay))) .* (x>=xresp) + (x<xresp)) * A .* ((M * (x-(xbegin-1))) + 1)';
    ft = fittype(eq, 'independent', 'x', 'dependent', 'y');

    excludedVals = false(1,length(trace));
    excludedVals(1:xbegin-1) = true;
    outliers = excludedata(xData,yData,'indices',excludedVals);
    opts = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', ...
        [0 -1 -inf 2 xbegin xresp], 'Upper', [2^16 20 0 30 xbegin xresp],...
        'Exclude', outliers, 'Display', 'notify');
    opts.StartPoint = [A B M tau_decay xbegin xresp];
    
    [fitresult, ~] = fit(xData, yData, ft, opts);

    % Placebo Fitting
    [xData, yData] = prepareCurveData(x(1:xresp-3), double(trace(1:xresp-3)));
    outliers = excludedata(xData,yData,'indices',excludedVals(1:xresp-3));
    coeffs = coeffvalues(fitresult);
    opts.Exclude = outliers;
    opts.Lower = [coeffs(1), -1, coeffs(3), 2,  coeffs(5), 10];
    opts.Upper = [coeffs(1), 20, coeffs(3), 30, coeffs(5), 10];
    opts.StartPoint = [coeffs(1) B coeffs(3) tau_decay coeffs(5) 10];
    [pseudoFit, ~] = fit(xData, yData, ft, opts);
end

function [fitresult, pseudoFit] = fitDoubleExponentialTrace(trace, xbegin, xresp)
    x = (1:length(trace));
    [xData, yData] = prepareCurveData(x, double(trace)); 

    A = median(trace(xbegin:xresp-3));
    B = mean(trace(xresp: xresp+2)) / mean(trace(xresp-4:xresp-2));
    tau_bleach = 10;
    tau_decay = 5;
    yoffset = 1;
    
    eq =  '((1 + (B * exp (-(x-xresp) ./ tau_decay))) .* (x>=xresp) + (x<xresp)) * A .* ((exp (-(x-(xbegin-1)) ./ tau_bleach)) .* (1-yoffset) + yoffset) ';
    ft = fittype(eq, 'independent', 'x', 'dependent', 'y');

    excludedVals = false(1,length(trace));
    excludedVals(1:xbegin-1) = true;
    outliers = excludedata(xData,yData,'indices',excludedVals);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', ...
        [0 -1 3 2 xbegin xresp 0], 'Upper',[2^16 20 50 30 xbegin xresp 1],...
        'Exclude', outliers, 'Display', 'notify');
    opts.StartPoint = [A B tau_bleach tau_decay xbegin xresp yoffset];
    
    [fitresult] = fit(xData, yData, ft, opts);

    % Placebo Fitting
    [xData, yData] = prepareCurveData(x(1:xresp-3), double(trace(1:xresp-3)));
    outliers = excludedata(xData,yData,'indices',excludedVals(1:xresp-3));
    coeffs = coeffvalues(fitresult);
    opts.Exclude = outliers;
    opts.Lower = [coeffs(1), -1, coeffs(3), 2, coeffs(5), 10, coeffs(7)];
    opts.Upper = [coeffs(1), 20, coeffs(3), 30, coeffs(5), 10, coeffs(7)];
    opts.StartPoint = [coeffs(1) B coeffs(3) tau_decay coeffs(5) 10 coeffs(7)];
    [pseudoFit, ~] = fit(xData, yData, ft, opts);
end

function normalizedTrace = normalizeTrace(trace, fitresult)
    coeffvals = coeffvalues(fitresult);
    normalizedTrace = trace / coeffvals(1);
end

function [correctedTrace, normalizedTrace] = getLineCorrectionAndNormalization(trace, fitresult, xbegin)
    x = (0:length(trace)-1);
    coeffvals = coeffvalues(fitresult);
    correctionVector = (coeffvals(3) * (x-(xbegin-1))) + 1;
    correctionVector(1:xbegin-1) = 1;
    correctedTrace = trace ./ correctionVector;
    normalizedTrace = correctedTrace / coeffvals(1);

%     figure('Position',[473,386,1009,592]);
%     sgtitle(["Weak Bleaching", "A = " + round(coeffvals(1),2) + ", B = " + round(coeffvals(2),2) + ", M = " + coeffvals(3) ...
%         + ", tau_decay = " + coeffvals(4)], 'Interpreter','none')
%     subplot(5,1,1); plot(trace); axis tight; grid minor; title("Origianl Trace")
% %     subplot(5,1,2); plot((xbegin+1:30), fitresult((xbegin+1:30)), 'r'); hold on; plot(1:30,trace); axis tight; grid minor; title("Trace + Fitting")
%     subplot(5,1,2); plot((1:30), fitresult((1:30)), 'r'); hold on; plot(1:30,trace); axis tight; grid minor; title("Trace + Fitting")
%     subplot(5,1,3); plot(correctionVector, 'k'); axis tight; grid minor; title("Correction Vector")
%     subplot(5,1,4); plot(correctedTrace); axis tight; grid minor; title("Corrected Trace")
%     subplot(5,1,5); plot(normalizedTrace); axis tight; grid minor; title("Normalized Trace")
end

function [correctedTrace, normalizedTrace] = getExponentialCorrectionAndNormalization(trace, fitresult, xbegin)
    x = (0:length(trace)-1);
    coeffvals = coeffvalues(fitresult);
    correctionVector = (exp (- (x-(xbegin-1)) / coeffvals(3))) * (1 - coeffvals(7)) + coeffvals(7);
    correctionVector(1:xbegin-1) = 1;
    correctedTrace = trace ./ correctionVector;
    normalizedTrace = correctedTrace / coeffvals(1);

%     figure('Position',[473,386,1009,592]);
%     sgtitle(["Strong Bleaching", "A = " + round(coeffvals(1),2) + ", B = " + round(coeffvals(2),2) + ", tau_bleach = " + coeffvals(3) ...
%         + ", tau_decay = " + coeffvals(4) + ", yoffset = " + round(coeffvals(7),2)], 'Interpreter','none')
%     subplot(5,1,1); plot(trace); axis tight; grid minor; title("Origianl Trace")
% %     subplot(5,1,2); plot((xbegin+1:30), fitresult((xbegin+1:30)), 'r'); hold on; plot(1:30,trace); axis tight; grid minor; title("Trace + Fitting")
%     subplot(5,1,2); plot((1:30), fitresult((1:30)), 'r'); hold on; plot(1:30,trace); axis tight; grid minor; title("Trace + Fitting")
%     subplot(5,1,3); plot(correctionVector, 'k'); axis tight; grid minor; title("Correction Vector")
%     subplot(5,1,4); plot(correctedTrace); axis tight; grid minor; title("Corrected Trace")
%     subplot(5,1,5); plot(normalizedTrace); axis tight; grid minor; title("Normalized Trace")
end

% Deprecated
function [correctedTrace, normalizedTrace, lineFit] = getSimpleLineCorrection(trace, xbegin)
    x = (1:length(trace));
    eq =  'A * ((M * x) + 1)';
    ft = fittype(eq, 'independent', 'x', 'dependent', 'y'); 

    [xData, yData] = prepareCurveData(x, double(trace)); 
    excludedVals = false(1,length(trace));
    excludedVals(1:xbegin) = true;
    outliers = excludedata(xData,yData,'indices',excludedVals);
    opts = fitoptions('Exclude', outliers);
%     opts.StartPoint = [mean(trace) 0];
    
    lineFit = fit(x',double(trace)',ft, opts);
    coeffvals = coeffvalues(lineFit);
    correctionVector = (coeffvals(2) * x) + 1;
    correctedTrace = trace ./ correctionVector;
    correctionVector = coeffvals(1) * ((coeffvals(2) * x) + 1);
    normalizedTrace = trace ./ correctionVector;
end

% Deprecated
function [correctedTrace, normalizedTrace, expFit] = getSimpleExpCorrection(trace, xbegin)
    x = (1:length(trace));
    eq =  'a * exp(-x * b)';
    ft = fittype(eq, 'independent', 'x', 'dependent', 'y'); 

    [xData, yData] = prepareCurveData(x, double(trace)); 
    excludedVals = false(1,length(trace));
    excludedVals(1:xbegin) = true;
    outliers = excludedata(xData,yData,'indices',excludedVals);
    opts = fitoptions('Method', 'NonlinearLeastSquares', 'Exclude', outliers, 'Lower', ...
        [0 0], 'Upper',[2^16 inf]);
    opts.StartPoint = [mean(trace) 0];

    expFit = fit(x',double(trace)',ft, opts);
    coeffvals = coeffvalues(expFit);
    correctionVector = exp(coeffvals(2) * x);
    correctedTrace = trace ./ correctionVector;
    correctionVector = coeffvals(1) * exp(coeffvals(2) * x);
    normalizedTrace = trace ./ correctionVector;
end

% TO CHECK: what to do with pseudofit
function traceCoeffs = getFittingCoefficients(fitresult, pseudoFit, baselineDecayPerc, fitType)
    coef = coeffvalues(fitresult);
    pseudoCoef = coeffvalues(pseudoFit);
    switch fitType
        case 'lineExp'
            traceCoeffs = [2, baselineDecayPerc, coef(1:2), coef(4), NaN(1,2), coef(3), pseudoCoef(2)];
        case 'DoubleExp'
            traceCoeffs = [3, baselineDecayPerc, coef(1:2), coef(4), coef(3), coef(7), NaN(1), pseudoCoef(2)];
        case 'expResponse'
            traceCoeffs = [1, baselineDecayPerc, coef(1:3), NaN(1,3), pseudoCoef(2)];
    end
end