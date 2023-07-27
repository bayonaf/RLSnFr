function signalImage = getSignalImage(deltaFMatrix, coordinate, submatSize, imageSize, triggerTime)
    
    [limitLeft, limitRight, limitUp, limitDown, isExceded] = getSubmatrixLimits(coordinate, fix(submatSize/2), imageSize);
    if isExceded == 1
        signalImage =...
            getPropperMatrix(deltaFMatrix(:,:,triggerTime),limitLeft, limitRight, limitUp, limitDown, submatSize, imageSize);
    else
        signalImage = deltaFMatrix(limitUp:limitDown, limitLeft:limitRight,triggerTime);
    end
end

function [limitLeft, limitRight, limitUp, limitDown, isExceded] = getSubmatrixLimits(coordinate, submatSize, imageSize)
    %GETSUBMATRIXLIMITS Summary of this function goes here
        isExceded = 0;
        limitLeft = coordinate(1, 2) - submatSize;
        if limitLeft < 1
            isExceded = 1;
        end
        limitRight = coordinate(1, 2) + submatSize;
        if limitRight > imageSize
            isExceded = 1;
        end
        limitUp = coordinate(1, 1) - submatSize;
        if limitUp < 1
            isExceded = 1;
        end
        limitDown = coordinate(1, 1) + submatSize;
        if limitDown > imageSize
            isExceded = 1;
        end
end

function submatBase = getPropperMatrix(deltaFMatrix, limitLeft, limitRight, limitUp, limitDown, sideSize, imageSize)
    
    [realLeft, newLeft, realRight, newRight, realUp, newUp, realDown, newDown] = ...
    deal(limitLeft, 1, limitRight, sideSize, limitUp, 1, limitDown, sideSize);
    submatBase = zeros(sideSize);
    if limitLeft < 1
        realLeft = 1;
        realRight = limitRight;
        newLeft =  2+(-limitLeft);
        newRight = sideSize;
    end
    if limitRight > imageSize
        realRight = imageSize;
        realLeft = limitLeft;
        newRight = sideSize-(limitRight - imageSize);
        newLeft = 1;
    end
    if limitUp < 1
        realUp = 1;
        realDown = limitDown;
        newUp = 2+(-limitUp);
        newDown = sideSize;
    end
    if limitDown > imageSize
        realDown = imageSize;
        realUp = limitUp;
        newDown = sideSize-(limitDown-imageSize);
        newUp = 1;
    end
     submatBase(newUp:newDown, newLeft:newRight) = deltaFMatrix(realUp:realDown, realLeft:realRight);
end





