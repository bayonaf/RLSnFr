function [defaultKernel, allDFSignalsCellVector, allCoordinatesCellVector,...
    normalizedCoeffs, averagedCoeffs, averagedCalibratedCoeffs, coeffsError, coeffsCV,...
    ampFunction, ampCalibratedFunction, sgmFunction, sgmClaibratedFunction, meanX,...
    meanCalibX, meanY, meanCalibY, gofAmp, gofSgm] = ...
    getCustomizedKernel(deltaFCellArray, coordinatesCellarray, sideSize, nrOfFrames,...
    kernelSideSieze, kernelNrOfFrames)
%   GETCUSTOMIZEDKERNEL: This function calculates a customized kernel for deconvolution, 
%   based on the previously known ROIs and signals of the stack 
%   CustomizedKernel = getCustomizedKernel(deltaFCellArray, coordinatesArray, sideSieze, nrOfFrames) 
%   returns a sideSieze x sideSieze x (2*nrOfFrames+1) matrix, with a
%   gaussian kernel centered in frame nrOfFrames+1
%   
%   Input:
%
%   IMAGESCELLARRAY:    Cell array containing all the 3D matrix stacks
%
%   DELTAFCELLARRAY:    Cell array containing all the 3D DF matrix stacks
%
%   COORDINATESARRAY:   Contains all the arrays with the corresponding
%   ROIs coordinates of every stack in DELTAFCELLARRAY
%
%   SIDESIZE:           Numerical value of the x,y sizes of the customized
%   kernel
%   
%   NROFFRAMES:         Numerical value of the number of frames that commposes 
%   the customized kernel gaussian. the kernel is centered in 3 dimensions, 
%   so is the number of frames before the gaussian starts in the kernel matrix   
%
%   Output:
%
%   KERNEL:             Customized kernel for deconvolution

%   1. All 3D submatriz containing signals defined by coordinates are
%   extracted
    disp("Calculating submatrix... ")
    [allSignalsCellVector, allCoordinatesCellVector] = deal(cell(length(deltaFCellArray),1));
    for i = 1: length(deltaFCellArray)
        if isequal(length(coordinatesCellarray), 1) % Only 1 coordinates file
            [allSignalsCellVector{i,1}, allCoordinatesCellVector{i,1}] = ...
                getAllSignalsSubMatrix(deltaFCellArray{i}, coordinatesCellarray{1,1}, sideSize, nrOfFrames);
        else % Many coordinates files
            [allSignalsCellVector{i,1}, allCoordinatesCellVector{i,1}] = ...
                getAllSignalsSubMatrix(deltaFCellArray{i}, coordinatesCellarray{i,1}, sideSize, nrOfFrames);
        end
    end
    
%   2. All 2D matrix containing the signals in time event zero are extracted 
%   to enhance propper localization in the ROI by gaussian fitting. This will
%   update the ROIs coordinates
    disp("Getting 2D signals in time event... ")
    allTimeEventZeroMatrix = cell(length(deltaFCellArray),1);
    for i = 1: length(deltaFCellArray)
        allTimeEventZeroMatrix{i,1} = getEventsInTimeEventFromSubmatrix(allSignalsCellVector{i,1}, allCoordinatesCellVector{i,1});
    end
    disp("Fitting signals with gaussian... ")
    allFitResults = cell(length(deltaFCellArray),1);
    for i = 1: length(deltaFCellArray)
        stackSize = size(deltaFCellArray{i,1}, 1);
        [allCoordinatesCellVector{i,1}, allFitResults{i,1}] = ...
            getCorrectedCoordinatesByGaussianFitting(allTimeEventZeroMatrix{i,1}, allCoordinatesCellVector{i,1}, stackSize);
    end
    
%   3. With the updated coordinates, the signals are extracted, and based 
%   on them, gaussian coefficients are calculated to produce a gaussian kernel
    disp("Getting Coefficients")
    allDFSignalsCellVector = cell(length(deltaFCellArray),1);
    for i = 1: length(deltaFCellArray)
        [allDFSignalsCellVector{i,1}, allCoordinatesCellVector{i,1}] = ...
            getAllSignalsSubMatrix(deltaFCellArray{i,1}, allCoordinatesCellVector{i,1}, sideSize, nrOfFrames);
    end
    [all2DGaussianFits, all2DGaussianCoefs] = getGaussianCoefficients(allDFSignalsCellVector, allCoordinatesCellVector);
    
%   4. Coefficients are normalized and averaged, and error in STD and Coefficients of variation are calculated    
%     disp("Averaging coefficients")
    [normalizedCoeffs, averagedCoeffs, averagedCalibratedCoeffs, coeffsError, coeffsCV] = ...
        getCoefficientsInformation(all2DGaussianCoefs);

%   5. Amplitude and sigma coefficients are fitted to find the formula they describe    
    disp("Getting new coefficients from fitting")
    [ampFunction, ampCalibratedFunction, sgmFunction, sgmClaibratedFunction,...
        meanX, meanCalibX, meanY, meanCalibY, gofAmp, gofSgm] = ...
        fitCoeficients(averagedCoeffs, averagedCalibratedCoeffs);

%   6. Knowing the formulas, a gaussian kernel is created
    disp("Generating kernel")
    defaultKernel = getSyntethicKernel(ampFunction, sgmFunction, kernelSideSieze, kernelNrOfFrames);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [signalsSubmatrixVector, coordinatesArray] = getAllSignalsSubMatrix(deltaFMatrix, coordinatesArray, sideSize, nrOfFrames)
%   GETALLSIGNALSSUBMATRIX Generates a cell array containing 3D submatrix of
%   every event defined by the coordinates of size submatSize x submatSize
%   x 2*nrOfFrames (the signals are centered in the x, y and frame coordinates of the submatrix. 
%   The fuction also update the coordinates vector with the new central
%   frame cooedinate where the signal starts

    if mod(sideSize, 2) == 0
        warning("Submatrix dimensions must be odd")
    end
    [~, imageSize, imageDeep] = size(deltaFMatrix);
    signalsSubmatrixVector = cell(size(coordinatesArray,1) ,1);
    f = waitbar(0,'Calculating Submatrix... 0%');
    for i = 1: size(coordinatesArray, 1)
        [downLimit, upLimit] = deal (coordinatesArray(i,3)-nrOfFrames,...
            coordinatesArray(i,3)+nrOfFrames); % Limit in frame axis where the signals will be substracted
        if nrOfFrames  < imageDeep % If the new limit exceed the limits of the original matrix, originals are keeped.
            if downLimit < 1
               if upLimit > imageDeep
                   [downLimit, upLimit] = deal (1, imageDeep);
               else
                   downLimit = 1;
                   coordinatesArray(i,4) = nrOfFrames-1;
               end
            else
                if upLimit > imageDeep
                   upLimit = imageDeep;
                   coordinatesArray(i,4) = nrOfFrames+1;
                else
                   coordinatesArray(i,4) = nrOfFrames+1;
                end
            end
        else
            [downLimit, upLimit] = deal (1, imageDeep);
        end
        coords = [coordinatesArray(i,2), coordinatesArray(i,1)];
        [limitLeft, limitRight, limitUp, limitDown, isExceded] =...
            getSubmatrixLimits(coords, fix(sideSize/2), imageSize); % the spatial limits where the signal values will be copied
        tempSubMat = zeros(limitRight-limitLeft+1, limitDown-limitUp+1, upLimit-downLimit+1);
        cont = 1;
        for j = downLimit: upLimit % Submatrix containing the signals are extracted
            if isExceded == 1
                tempSubMat(:,:,cont) =...
                    getPropperMatrix(deltaFMatrix(:,:,j),limitLeft, limitRight, limitUp, limitDown, sideSize, imageSize);
            else
                tempSubMat(:,:,cont) =...
                    deltaFMatrix(limitUp:limitDown, limitLeft:limitRight,j);
            end
            cont = cont + 1;
        end
        signalsSubmatrixVector{i} = tempSubMat; % All the signal submatrix
        waitbar(i/size(coordinatesArray, 1), f, 'Calculating Submatrix... ('...
            +string(i*100/size(coordinatesArray, 1))+'%)');
    end
    close(f)
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

function [limitLeft, limitRight, limitUp, limitDown, isExceded] = getSubmatrixLimits(coordinate, submatSize, imageSize)
    isExceded = 0;
    limitLeft = coordinate(1, 2) - submatSize;
    if limitLeft < 1
        warning("Left boundaries exceeded. Submatrix will contain the possible values of the images and values out of bounds will be 0")
        isExceded = 1;
    end
    limitRight = coordinate(1, 2) + submatSize;
    if limitRight > imageSize
        warning("Right boundaries exceeded. Submatrix will contain the possible values of the images and values out of bounds will be 0")
        isExceded = 1;
    end
    limitUp = coordinate(1, 1) - submatSize;
    if limitUp < 1
        warning("Upper boundaries exceeded. Submatrix will contain the possible values of the images and values out of bounds will be 0")
        isExceded = 1;
    end
    limitDown = coordinate(1, 1) + submatSize;
    if limitDown > imageSize
        warning("Bottom boundaries exceeded. Submatrix will contain the possible values of the images and values out of bounds will be 0")
        isExceded = 1;
    end
end

function signalsInTimeEventArray = getEventsInTimeEventFromSubmatrix(signalsSubmatrixVector, coordinatesArray)
%   GETEVENTSINTIMEZEROFROMCLEANSUBMATRIX Returns the 2D matrix in time event
%   zero of every signal submatrix
    signalsInTimeEventArray = cell(length(signalsSubmatrixVector), 1);
    for i = 1: size(signalsSubmatrixVector, 1)
        timeCoordinate = coordinatesArray(i,4);
        signalsInTimeEventArray{i, 1} = signalsSubmatrixVector{i,1}(:,:,timeCoordinate);
    end
end

function [coordinatesVector, fitResults] = ...
    getCorrectedCoordinatesByGaussianFitting(signalsSubmatrixVector, coordinatesVector, stackSize)
%GETCORRECTEDCOORDINATESBYGAUSSIANFITTING Summary of this function goes here
    fitResults = cell(length(signalsSubmatrixVector), 1);
    f = waitbar(0,'Fitting signals with gaussian... 0%');
    for i = 1: length(signalsSubmatrixVector)
        signalsSubmatrixVector{i,1}(signalsSubmatrixVector{i,1} < 0) = 0; % For gaussian fitting all values must be over zero
        [fitResult, ~] = createFit(signalsSubmatrixVector{i,1}); % Gaussian Fitting executed
        
        centerPos = fix(size(signalsSubmatrixVector{i,1}, 1)/2)+1;
        if (fitResult.x0 > 1) && (fitResult.x0 < stackSize)
            if fitResult.x0 > centerPos
                coordinatesVector(i, 1) = coordinatesVector(i, 1) + round(fitResult.x0) - centerPos;
            end
            if fitResult.x0 < centerPos
                coordinatesVector(i, 1) = coordinatesVector(i, 1) - (centerPos - round(fitResult.x0));
            end
        end
        if (fitResult.y0 > 1) && (fitResult.y0 < stackSize)
            if fitResult.y0 > centerPos
                coordinatesVector(i, 2) = coordinatesVector(i, 2) + round(fitResult.y0) - centerPos;
            end
            if fitResult.y0 < centerPos
                coordinatesVector(i, 2) = coordinatesVector(i, 2) - (centerPos - round(fitResult.y0));
            end
        end
        fitResults{i} = fitResult;
        waitbar(i/length(signalsSubmatrixVector), f,'Fitting signals with gaussian... ('+string(i*100/length(signalsSubmatrixVector))+'%)');
    end
    close(f)
end

function [fitResult, gof] = createFit(submatrix)
    [x, y] = deal((1:size(submatrix, 1)), (1:size(submatrix, 2)));
    [xData, yData, zData] = prepareSurfaceData(x, y, submatrix);

    % Set up fittype and options.
    ft = fittype( 'Amp*exp(-(((x-x0)^2)/(2*sgmX^2) + ((y-y0)^2)/(2*sgmY^2)))', 'independent', {'x', 'y'}, 'dependent', 'z' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.MaxFunEvals = 5000;
    opts.MaxIter = 5000;
    opts.StartPoint = [100 5 5 10 10]; % Hay que calcular esto antes
    opts.TolFun = 1e-12;
    opts.TolX = 1e-12;

    % Fit model to data.
    [fitResult, gof] = fit([xData, yData], zData, ft, opts);
%     figure( 'Name', 'untitled fit 1' );
%     h = plot( fitResult, [xData, yData], zData );
%     legend( h, 'untitled fit 1', 'M vs. x, x', 'Location', 'NorthEast', 'Interpreter', 'none' );
%     % Label axes
%     xlabel( 'x', 'Interpreter', 'none' );
%     ylabel( 'x', 'Interpreter', 'none' );
%     zlabel( 'M', 'Interpreter', 'none' );
%     grid on
%     view( -89.9, 1.6 );
end

function [allGaussianFits, allGaussianCoefs] = getGaussianCoefficients(stacksSignals, allCoordinatesCellVector)
% GET2DGAUSSIANCOEFFICIENTS Summary of this function goes here
    disp("Calculating Gaussian Coefficients for kernel")
    [allGaussianFits, allGaussianCoefs] = deal(cell(0, 1));
%     f = waitbar(0,"Fitting signals with gaussian... 0%");
    for i = 1: length(stacksSignals) % Number of stacks
        [gaussianFits, gaussianCoeffs] = deal(cell(0, 1));
        for j = 1: length(stacksSignals{i,1}) % Number of signals per stack
            numFrames = size(stacksSignals{i,1}{j,1}, 3);
            if isequal(length(allCoordinatesCellVector), 1)
            	timeEvent = allCoordinatesCellVector{1,1}(j,4);
            else
                timeEvent = allCoordinatesCellVector{i,1}(j,4);
            end
            fits = cell(0, 1);
            coeffs = zeros(numFrames-3 - timeEvent, 5);
            for k = timeEvent: numFrames-3 % 2D matrix after the synapse release
                [fitresult, gof] = createFit(stacksSignals{i,1}{j,1}(:,:,k));
                f = struct("FitResult", fitresult, "gof", gof);
                fits{k-timeEvent+1,1} = f;
                coeffs(k-timeEvent+1,:) = coeffvalues(fitresult);
%                 waitbar(i/length(numFrames), f,"Fitting signals with gaussian for Stack " + i + ", signal ");
            end
            gaussianFits{j} = fits;
            gaussianCoeffs{j} = coeffs;
        end
        allGaussianFits{i} = gaussianFits;
        allGaussianCoefs{i} = gaussianCoeffs;
    end
%     close(f)
    
end

function [normalizedCalibratedCoeffs, averagedCoeffs, averagedCalibratedCoeffs, coeffsError, coeffsCV] = getCoefficientsInformation(coeeficients)
% GETCOEFFICIENTSINFORMATION Summary of this function goes here
% coefficients with correct calibration in µm
    calibratedCoeffs = coeeficients;
    for i = 1: length(calibratedCoeffs)
        for j = 1: length(calibratedCoeffs{i})
            calibratedCoeffs{i}{j}(:,2:4) = calibratedCoeffs{i}{j}(:,2:4) * 0.254;
        end
    end

    % Coefficients Normalization
    normalizedCalibratedCoeffs = cell(0,1);
    normalizedCoeffs = cell(0,1);
    cont = 1;
    for i = 1: length(calibratedCoeffs) % Iterating over stacks
        for j = 1: length(calibratedCoeffs{i}) % Iterating over signals per stack
            
            temp = zeros(length(coeeficients{i}{j}), 4);
            temp(:,1) = coeeficients{i}{j}(:,1) ./ coeeficients{i}{j}(1,1);
            temp(:,2:5) = coeeficients{i}{j}(:,2:5);
            normalizedCoeffs{cont, 1} = temp;
            
            temp = zeros(length(calibratedCoeffs{i}{j}), 4);
            temp(:,1) = calibratedCoeffs{i}{j}(:,1) ./ calibratedCoeffs{i}{j}(1,1);
            temp(:,2:5) = calibratedCoeffs{i}{j}(:,2:5);
            normalizedCalibratedCoeffs{cont, 1} = temp;
            cont = cont + 1;
        end
    end
    normalizedCoeffs = cat(3, normalizedCoeffs{:});
    normalizedCalibratedCoeffs = cat(3, normalizedCalibratedCoeffs{:});
    
    % Coefficients Average
    averagedCoeffs = mean(normalizedCoeffs, 3);
    averagedCalibratedCoeffs =  mean(normalizedCalibratedCoeffs, 3);
    
    % Coefficients Error
    coeffsError = std(normalizedCalibratedCoeffs, [], 3)/sqrt(size(normalizedCalibratedCoeffs, 3));
    
    % Coefficients of Variation    
    coeffsCV = std(normalizedCalibratedCoeffs, [], 3)./averagedCalibratedCoeffs;    
end

function [ampFunction, ampCalibratedFunction, sgmFunction, sgmClaibratedFunction, meanX, meanCalibX, meanY, meanCalibY, gofAmp, gofSgm]...
    = fitCoeficients(averagedCoeffs, averagedCalibratedCoeffs)
    X = (10: 20: 130);
    [ampFunction, ~] = createInverseFit(X', averagedCoeffs(:,1)');
    [sgmFunction, ~] = createPoly2Fit(X', averagedCoeffs(:,2)');
    [ampCalibratedFunction, gofAmp] = createInverseFit(X', averagedCalibratedCoeffs(:,1)');
    [sgmClaibratedFunction, gofSgm] = createPoly2Fit(X', averagedCalibratedCoeffs(:,2)');
    meanX = mean(averagedCoeffs(:,3));
    meanY = mean(averagedCoeffs(:,4));
    meanCalibX = mean(averagedCalibratedCoeffs(:,3));
    meanCalibY = mean(averagedCalibratedCoeffs(:,4));
end

function [fitresult, gof] = createInverseFit(X, A)
    [xData, yData] = prepareCurveData( X, A );

    % Set up fittype and options.
    ft = fittype( 'a1 * [exp(-x/tauOff) - exp(-x/tauOn)]', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [1 100 10];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
end

function [fitresult, gof] = createPoly2Fit(X, A)
    [xData, yData] = prepareCurveData( X, A );

    % Set up fittype and options.
    ft = fittype( 'poly1' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
end


function kernel = getSyntethicKernel(ampFunction, sgmFunction, kernelSideSieze, kernelNrOfFrames)
    % GETSYNTETHICKERNEL Summary of this function goes here
    framesNumber = (kernelNrOfFrames*2+1); % Kernel size in time
    kernelNrOfFrames = fix(framesNumber/2)+1; % Kernel center in time;
    center = fix(kernelSideSieze/2)+1;
    kernel = zeros(kernelSideSieze, kernelSideSieze, framesNumber); % Kernel Base
    numberOfSignals = framesNumber - kernelNrOfFrames + 1; % Number of gaussian curves must be generated
    X = (10: 20: framesNumber*10); % Evaluation interval for the functions
    
    % Calculating coefficients from regression functions
    ampCoeffs = ampFunction(X);
    sgmCoeffs = sgmFunction(X);
    
    % Filling kernel's base with 2D gausssians according to coefficients
    for i = 1: numberOfSignals
%         kernel(:,:,centralFrame) = gauss2d(ampCoeffs(i), meanX, meanY, sgmCoeffs(i), sgmCoeffs(i), kernelSize);
        kernel(:,:,kernelNrOfFrames) = gauss2d(ampCoeffs(i), center, center, sgmCoeffs(i), sgmCoeffs(i), kernelSideSieze);
        kernelNrOfFrames = kernelNrOfFrames + 1;
    end
end

function g = gauss2d(amp, x0, y0, sigmax, sigmay, tam)
    [x, y] = meshgrid(1:1:tam, 1:1:tam);
    xExpression = ((x - x0) .^2) / (2.0 * (sigmax .^2));
    yExpression = ((y - y0) .^2) / (2.0 * (sigmay .^2));
    g = amp * exp(-((xExpression) + (yExpression)));
end