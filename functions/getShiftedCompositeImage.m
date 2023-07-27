function RGB = getShiftedCompositeImage(imagesCellArray)
%GETSHIFTEDCOMPOSITEIMAGE Summary of this function goes here
    stacksNumber = length(imagesCellArray);
    cmap = hsv(stacksNumber) / (stacksNumber/2);
    [m,n,~] = size(imagesCellArray{1});
    
    allFrames = [];
    for i = 1 : stacksNumber
        allFrames =  cat(3, allFrames, mean(imagesCellArray{i}(:,:,4:17), 3));
    end

    RGB = zeros(m, n, 3);
    for i = 1 : m
        for j = 1 : n
            [sR, sG, sB] = deal(0);
            for f = 1 : stacksNumber
                sR = sR + allFrames(i,j,f)/stacksNumber * cmap(f,1);
                sG = sG + allFrames(i,j,f)/stacksNumber * cmap(f,2);
                sB = sB + allFrames(i,j,f)/stacksNumber * cmap(f,3);
            end
            RGB(i,j,1) = sR;
            RGB(i,j,2) = sG;
            RGB(i,j,3) = sB;
        end
    end
end

function [RGB, shifting] = correctImageBackground(RGB, regCoeffs)
    [avgR, avgG, avgB] = deal(mean(RGB(:,:,1), 'all'), mean(RGB(:,:,2), 'all'), mean(RGB(:,:,3), 'all'));
    [shifting.MXn, shifting.MXp, shifting.MYn, shifting.MYp] = deal(0);
    
    % Negative values in X
    N = regCoeffs(2,:);
    idx = N < 0;
    neg = N(idx);
    % Positive values in X
    pos = N(~idx);

    if ~isempty(neg)
        MX = ceil(sqrt(max(neg.^2)));
        RGB(end-MX:end, :, 1) = avgR;
        RGB(end-MX:end, :, 2) = avgG;
        RGB(end-MX:end, :, 3) = avgB;
        shifting.MXn = MX;
    end
    if ~isempty(pos)
        MX = ceil(max(pos));
        RGB(1:MX, :, 1) = avgR;
        RGB(1:MX, :, 2) = avgG;
        RGB(1:MX, :, 3) = avgB;
        shifting.MXp = MX;
    end

    % Negative values in Y
    N = regCoeffs(1,:);
    idx = N < 0;
    neg = N(idx);
    % Positive values in Y
    pos = N(~idx);

    if ~isempty(neg)
        MY = ceil(sqrt(max(neg.^2)));
        RGB(:, end-MY:end, 1) = avgR;
        RGB(:, end-MY:end, 2) = avgG;
        RGB(:, end-MY:end, 3) = avgB;
        shifting.MYn = MY;
    end
    if ~isempty(pos)
        MY = ceil(max(pos));
        RGB(:, 1:MY, 1) = avgR;
        RGB(:, 1:MY, 2) = avgG;
        RGB(:, 1:MY, 3) = avgB;
        shifting.MYp = MY;
    end
end
