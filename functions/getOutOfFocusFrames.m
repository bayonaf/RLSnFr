function res = getOutOfFocusFrames(imageStack, stdFactor)
%GETOUTOFFOCUSFRAMES Summary of this function goes here
    
    SQRTSSE = zeros(1, size(imageStack, 3) - 1);
    f = waitbar(0,"Please Wait");
    tic
    steps = size(imageStack, 3) - 1;
    for i = 1 : steps
        SQRTSSE(i) = norm(imageStack(:,:,i) - imageStack(:,:,i+1)); % sqrt(sse(A-B))
        waitbar(i / steps , f, string(round(i * 100 / steps, 0)) + "%");
    end
    meanStack = reshape(mean(imageStack, [2 1], 'omitnan'), [1,size(imageStack, 3)]);
    SQRTSSE = SQRTSSE ./ meanStack(2:end);
    delete(f)
    toc
    
    SQRTSSEnormalized = SQRTSSE - mean(SQRTSSE);
    STD = std(SQRTSSEnormalized);

    [pksPos, locsPos] = findpeaks(SQRTSSEnormalized, "MinPeakHeight", stdFactor * STD);
    [pksNeg, locsNeg] = findpeaks(-SQRTSSEnormalized, "MinPeakHeight", stdFactor * STD);

    res.SQRTSSE = SQRTSSE;
    res.meanStack = meanStack;
    res.SQRTSSEnormalized = SQRTSSEnormalized;
    res.STD = STD;
    res.pksPos = pksPos;
    res.locsPos = locsPos;
    res.pksNeg = pksNeg;
    res.locsNeg = locsNeg;
end

