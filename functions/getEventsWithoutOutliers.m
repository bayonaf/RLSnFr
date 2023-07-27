function correctedEventPoints = getEventsWithoutOutliers(experimentType, eventsPoints, regCoeffs, maxSideSize, windowOffset, windowLength)
%EXCLUDEOULIEREVENTS Summary of this function goes here
    
    % Get maximum Shiftings
    displacement = getMaxShifting(regCoeffs);
    % Delete the outlier points
    correctedEventPoints = correctevents(experimentType, eventsPoints, displacement, maxSideSize, windowOffset, windowLength);

end

function correctedEventPoints = correctevents(experimentType, eventsPoints, displacement, maxSideSize, windowOffset, windowLength)
    stacksNumber = length(eventsPoints);
    correctedEventPoints = cell(1);
    
    for i = 1: stacksNumber
        trialEvents = eventsPoints{i};

        % Delete events in Right Border
        X = trialEvents(:,2);
        idx = X >= maxSideSize - displacement.MaxRigthX;
        trialEvents(idx, :) = [];

        % Delete events in Left Border
        X = trialEvents(:,2);
        idx = X <= displacement.MaxLeftX;
        trialEvents(idx, :) = [];

        % Delete events in Bottom  Border
        Y = trialEvents(:,3);
        idx = Y >= maxSideSize - displacement.MaxBottomY;
        trialEvents(idx, :) = [];

        % Delete events in Left Border
        Y = trialEvents(:,3);
        idx = Y <= displacement.MaxTopY;
        trialEvents(idx, :) = [];

        correctedEventPoints{i} = trialEvents;
    end

    if experimentType == 1
        T = trialEvents(:,4);
        idx = T <= (2 * windowOffset) + windowLength;
        correctedEventPoints{1}(idx, :) = [];
    end
end

function displacement = getMaxShifting(regCoeffs)
    [displacement.MaxRigthX, displacement.MaxLeftX, displacement.MaxBottomY, ...
        displacement.MaxTopY] = deal(0);
    bias = 2;

    % Negative values in X
    N = regCoeffs(1,:);
    idx = N < 0;
    neg = N(idx);
    % Positive values in X
    pos = N(~idx);

    if ~isempty(neg)
        displacement.MaxRigthX = ceil(sqrt(max(neg.^2))) + bias;
    end
    if ~isempty(pos)
        displacement.MaxLeftX = ceil(max(pos)) + bias;
    end

    % Negative values in Y
    N = regCoeffs(2,:);
    idx = N < 0;
    neg = N(idx);
    % Positive values in Y
    pos = N(~idx);

    if ~isempty(neg)
        displacement.MaxBottomY = ceil(sqrt(max(neg.^2))) + bias;
    end
    if ~isempty(pos)
        displacement.MaxTopY = ceil(max(pos)) + bias;
    end
end

