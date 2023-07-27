function baselinesAVG = getSigmoidMask(baselinesAVG, thresh, maxBranchesDistNm, pixelNanoEquiv, sigmoidVector)
%GETSIGMOIDMASK Summary of this function goes here
    
    pixelNm = maxBranchesDistNm / pixelNanoEquiv;
    baselinesAVG(baselinesAVG < thresh) = 0;
    baselinesAVG(baselinesAVG >= thresh) = 1;
    
    [row, col] = find(baselinesAVG == 0);
    pointsToCheck = [row, col];
    [row, col] = find(baselinesAVG == 1);
    branchesPoints = [row, col];
    
    [~, D] = knnsearch(branchesPoints, pointsToCheck);
    
    for i = 1 : length(pointsToCheck)
        if D(i) <= pixelNm % Less than 1000nm
            baselinesAVG(pointsToCheck(i,1), pointsToCheck(i,2)) = 0.5;
        end
    end
    
    [row, col] = find(baselinesAVG == 0.5);
    pointsToCheck = [row, col];
    
    [~, D] = knnsearch(branchesPoints, pointsToCheck);
    
    distancesVector = linspace(0, pixelNm, maxBranchesDistNm);
    
    closestValue = zeros(1,length(D));
    for i = 1 : length(D)
        closestIndex = length(find(distancesVector <= D(i)));
        closestValue(i) = sigmoidVector(closestIndex);
        baselinesAVG(pointsToCheck(i,1), pointsToCheck(i,2)) = closestValue(i);
    end
end

