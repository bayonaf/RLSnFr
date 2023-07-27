function [unbleachedImagesCellArray, meanPctDiff, ...
    allAvgBefore, allAvgAfterItrial, allAvgAfterWhole, allFits, allCoeffs]...
    = getImagesWithoutBleaching(imagesCellArray, xbegin, xresp)

% GETIMAGESWITHOUTBLEACHING This function calculates the exponential fitting
%   of the averaged brigthness of each image stack, normalize it, and
%   correct the bleaching
%   CLEANIMAGESCELLARRAY = getImagesWithoutBleaching(clearImagesCellArray)
%   return the modified cell array with the unbleached images
    xresp = xresp-1;
    stacksNumber = length(imagesCellArray);
    sSize = size(imagesCellArray{1},3);
    x = (0:sSize-1);
    
    % 1. Preallocating memory for results 
    allAvgBefore = zeros(stacksNumber, sSize); % Trial's frames averages before unbleaching
    allAvgAfterItrial = zeros(stacksNumber, sSize); % Trial's frames averages after unbleaching
    allAvgAfterWhole = zeros(stacksNumber, sSize); % Across trials averages after unbleaching
    allFits = zeros(stacksNumber, sSize);
    allCoeffs = zeros(stacksNumber, 6);
    unbleachedImagesCellArray = cell(1);
    
    % 2. Calculating initial points
    baseCoeffs = getBaseCoefficients(imagesCellArray, xbegin, xresp, stacksNumber);
    A = baseCoeffs(1);
    B = baseCoeffs(2);
    tau_bleach = baseCoeffs(3);
    tau_decay = baseCoeffs(4);
    yoffset = baseCoeffs(6);
    
    for i = 1 : stacksNumber
        % 3. Calculate means vectors and exp fitting
        meanStackBeforeCorrection = reshape(mean(imagesCellArray{i}, [2 1], 'omitnan'), [1,sSize]);
        
        [xData, yData] = prepareCurveData(x, double(meanStackBeforeCorrection)); % Data to fit
        
%         eq =  '(exp (-(x) / tau_bleach) * (1-yoffset) + yoffset) * (A * (((B * exp (-(x-xresp) / tau_decay)) + 1) * (x>=xresp) + (x<xresp)))';
        eq = '(exp (-(x) / tau_bleach) * (1-yoffset) + yoffset) * (A * (1 + (B * exp (- (x-xresp)/tau_decay) * (x>=xresp))))';
        ft = fittype(eq, 'independent', 'x', 'dependent', 'y');
        excludedVals = false(1,length(yData)); % values to be ignored in fitting
        excludedVals(1:xbegin) = true;
        outliers = excludedata(xData,yData,'indices',excludedVals);
        opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', ...
            [0 0 tau_bleach tau_decay xresp 0], 'Upper',[inf inf tau_bleach tau_decay xresp 1],...
            'MaxIter', 30000, 'Exclude', outliers, 'Display', 'Off');
        opts.StartPoint = [A B tau_bleach tau_decay xresp yoffset];
        [fitresult, ~] = fit(xData, yData, ft, opts);
        
        % 3. Calculating normalization factors vector from fitting results
        xData = (0:sSize-1);
        coeffvals = coeffvalues(fitresult);
        fitMeanStack = ((exp (- (xData) / coeffvals(3)) * (1 - coeffvals(6))) + coeffvals(6));
        fitVector = fitresult(xData);
        fitVector(1:xbegin) = NaN;
        allFits(i,:) = fitVector;
        allCoeffs(i,:) = coeffvals;
        
        % 4. Normalizing frames in stack according to normalization vector
        % (bleaching correction)
        for j = 1: sSize
            unbleachedImagesCellArray{i}(:,:,j) = single(imagesCellArray{i}(:,:,j)) / fitMeanStack(j);
        end
        
        % 5. Calculating corrected averages vector
        meanStackAfterCorrection = reshape(mean(unbleachedImagesCellArray{i}, [2 1], 'omitnan'), [1,sSize]);
      
        % Accumulating all stacks data
        allAvgBefore(i,:) = meanStackBeforeCorrection;
        allAvgAfterItrial(i,:) = meanStackAfterCorrection;
%         allROIInTrial(i,:) = reshape(unbleachedImagesCellArray{i}(224, 265, :),1,[]);
    end
    
    % 6. Across trials correction
    firstStackMean = mean(allAvgAfterItrial(1,:));
    for i = 1 : stacksNumber
        stackMean = mean(allAvgAfterItrial(i,:));
        for j = 1: sSize
            unbleachedImagesCellArray{i}(:,:,j) = ...
                unbleachedImagesCellArray{i}(:,:,j) / (stackMean / firstStackMean);
        end
        meanStackAfterWholeCorrection = reshape(mean(unbleachedImagesCellArray{i}, [2 1], 'omitnan'), [1,sSize]);
        allAvgAfterWhole(i,:) = meanStackAfterWholeCorrection;
    end

    % 7. Exponential fitting across trials
    AllStacksMeansBefore = mean(allAvgBefore,2);
    AllStacksMeansAfter = mean(allAvgAfterWhole,2);
    deltaSignal = abs(AllStacksMeansAfter - AllStacksMeansBefore);
    percentageDifference = deltaSignal ./ AllStacksMeansAfter; % Percent by element.
    meanPctDiff = mean(percentageDifference) * 100; % Average percentage over all elements.   
    
end

function baseCoeffs = getBaseCoefficients(imagesCellArray, xbegin, xresp, stacksNumber)
    
    sSize = size(imagesCellArray{1},3);
    x = (0:sSize-1);
    allAvgBefore = zeros(stacksNumber, sSize);
    for t = 1 : stacksNumber
        allAvgBefore(t,:) = reshape(mean(imagesCellArray{t}, [2 1], 'omitnan'), [1,sSize]);
    end
    allAveraged = mean(allAvgBefore,1);
    
    [xData, yData] = prepareCurveData(x,allAveraged); % Data to fit
    yoffset = 0.5;
    A = yData(xbegin);
    B = (yData(xresp+1) - yData(xresp))  / yData(xresp+1);
    if B < 0.003
        B = 0.003;
    end
    tau_bleach = -(xresp-1) / (log(yData(xresp-1) / yData(xbegin) - yoffset));
    tau_decay = 10;
    
%     eq = '(exp (-(x) / tau_bleach) * (1-yoffset) + yoffset) * (A * (((B * exp (-(x-xresp) / tau_decay)) + 1) * (x>=xresp) + (x<xresp)))';
    eq = '(exp (-(x) / tau_bleach) * (1-yoffset) + yoffset) * (A * (1 + (B * exp (- (x-xresp)/tau_decay) * (x>=xresp))))';
    ft = fittype(eq, 'independent', 'x', 'dependent', 'y');
    excludedVals = false(1,length(yData)); % values to be ignored in fitting
    excludedVals(1:xbegin) = true;
    outliers = excludedata(xData,yData,'indices',excludedVals);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', ...
        [0 0 0 0 xresp 0], 'Upper',[inf inf inf inf xresp 1],...
        'MaxIter', 30000, 'Exclude', outliers, 'Display', 'Off');
    opts.StartPoint = [A B tau_bleach tau_decay xresp yoffset];
    [fitresult, ~] = fit(xData, yData, ft, opts);
    baseCoeffs = coeffvalues(fitresult);
end
