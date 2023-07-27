function correctedEventPoints = getEventsWithoutOutliersFromRegistered(eventsPoints, registeredStacks, maxSideSize)
%GETEVENTSWITHOUTOUTLIERSFROMREGISTERED Summary of this function goes here

    % Calculate max. num. pixels moved
    bias = 2;
    [allRowZeros, allColZeros] = calculateMaxMoved(registeredStacks);
    allRowZeros = allRowZeros + bias;
    allColZeros = allColZeros + bias;

    stacksNumber = length(eventsPoints);
    correctedEventPoints = cell(1);
    
    for i = 1: stacksNumber
        trialEvents = eventsPoints{i};

        % Delete events in Right Border
        X = trialEvents(:,2);
        idx = X >= maxSideSize - allColZeros;
        trialEvents(idx, :) = [];

        % Delete events in Left Border
        X = trialEvents(:,2);
        idx = X <= allColZeros;
        trialEvents(idx, :) = [];

        % Delete events in Bottom  Border
        Y = trialEvents(:,3);
        idx = Y >= maxSideSize - allRowZeros;
        trialEvents(idx, :) = [];

        % Delete events in Left Border
        Y = trialEvents(:,3);
        idx = Y <= allRowZeros;
        trialEvents(idx, :) = [];

        correctedEventPoints{i} = trialEvents;
    end
end

function [allRowZeros, allColZeros] = calculateMaxMoved(registeredStacks)
    [allRowZeros, allColZeros] = deal(0);
    for i = 1 : length(registeredStacks)
        rowZeros = 0;
        for j = 1 : size(registeredStacks{i},3)
            temp = numel(find(all(registeredStacks{i}(:,:,j) == 0,1)));
            if temp > rowZeros
                rowZeros = temp;
            end
        end
        if rowZeros > allRowZeros
            allRowZeros = rowZeros;
        end
    
        colZeros = 0;
        for j = 1 : size(registeredStacks{i},3)
            temp = numel(find(all(registeredStacks{i}(:,:,j) == 0,2)));
            if temp > colZeros
                colZeros = temp;
            end
        end
        if colZeros > allColZeros
            allColZeros = colZeros;
        end
    end
end

