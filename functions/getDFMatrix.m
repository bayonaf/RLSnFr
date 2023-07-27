function varargout = getDFMatrix(varargin)
%   GETDFMatrix: returns 3D delta F stack and 2D background
%   [DELTAFMATRIX, background] = getDFMatrix(IMAGEMATRIX, downLimit, upLimit)
%   return the delta F of a given stack and the corresponding background
%   substracted. Background is an averaged matrix of the stack frames between 
%   downLimit and uplimit
%
%   IMAGEMATRIX:   3D numeric array
%
%   DOWNLIMIT:     Number of the starting frame that will be averaged to
%   calculate the background
%
%   UPLIMIT:       Number of the ending frame that will be averaged to
%   calculate the background

    switch varargin{1}
        case 0 % Triggered
            frames = size(varargin{2}, 3);
            A = varargin{2}(:,:,varargin{3} : varargin{4});
            baseline = mean(A, 3);
            DELTAFMATRIX = zeros(size(varargin{2}), 'single');
            DELTAFOVERF = zeros(size(varargin{2}), 'single');
            for i = 1: frames
                DELTAFMATRIX(:, :, i) = single(varargin{2}(:, :, i)) - single(baseline);
                DELTAFOVERF(:, :, i) = DELTAFMATRIX(:, :, i) ./ single(baseline);
            end
            DELTAFOVERF(isnan(DELTAFOVERF)) = 0;
            DELTAFOVERF(isinf(DELTAFOVERF)) = 0;

            varargout{1} = DELTAFMATRIX;
            varargout{2} = DELTAFOVERF;
            varargout{3} = baseline;
        case 1 % Espontaneous
            registeredImages = cat(3, varargin{2}{:});
            startFrame = varargin{4} + varargin{3};
            
            [DFImage, DFFImage] = deal(zeros(size(registeredImages)));
            sSize = size(registeredImages, 3);
            startPt = 1;
            
            for f = startFrame : sSize
                baseline = registeredImages(:,:, startPt:startPt + varargin{4});
                baselineMean = mean(baseline,3);
                DFImage(:,:,f)  = registeredImages(:,:,f) - baselineMean;
                DFFImage(:,:,f) = (registeredImages(:,:,f) - baselineMean) ./ baselineMean; 
                startPt = startPt + 1;
            end
            DFImage(:,:,1:startFrame)  = 0;
            DFFImage(:,:,1:startFrame) = 0;

            varargout{1} = DFImage;
            varargout{2} = DFFImage;
    end
end

function baseline = correctBaselineZeros(baseline)
    [row, col] = find(~baseline);
    for i = 1 : length(row)
        [limitLeft, limitRight, limitUp, limitDown, isExceded] = getSubmatrixLimits([row(i), col(i)], size(baseline,1));
        if isExceded == 0
            submat = baseline(limitUp:limitDown, limitLeft:limitRight);
        else
            submat = getPropperMatrix(baseline , limitLeft, limitRight, limitUp, limitDown, 3, size(baseline,1));
        end
        submat = submat(submat > 0);
        baseline(row, col) = mean(submat);
    end
end

function [limitLeft, limitRight, limitUp, limitDown, isExceded] = getSubmatrixLimits(coordinate, imageSize)
    isExceded = 0;
    limitLeft = coordinate(1, 2) - 1;
    if limitLeft < 1
        isExceded = 1;
    end
    limitRight = coordinate(1, 2) + 1;
    if limitRight > imageSize
        isExceded = 1;
    end
    limitUp = coordinate(1, 1) - 1;
    if limitUp < 1
        isExceded = 1;
    end
    limitDown = coordinate(1, 1) + 1;
    if limitDown > imageSize
        isExceded = 1;
    end
end

function submatBase = getPropperMatrix(matrixToInsert, limitLeft, limitRight, limitUp, limitDown, sideSize, imageSize)
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
    submatBase(newUp:newDown, newLeft:newRight) = matrixToInsert(realUp:realDown, realLeft:realRight);
end