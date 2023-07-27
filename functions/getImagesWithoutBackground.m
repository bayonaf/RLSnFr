function [clearImagesCellArray, avgs, pos] = getImagesWithoutBackground(imagesCellArray)
% GETIMAGESWITHOUTBACKGROUND This function returns the images without the
%   background noise
%   CLEARIMAGESCELLARRAY = getImagesWithoutBackground(imagesCellArray) 
%   return a cell array containing the images stacks without the average
%   value of the background of an image region selected manually by the
%   user with the mouse.
%
%   IMAGESCELLARRAY:   Cell array of 3D stacks

    stacksNumber = length(imagesCellArray);
    clearImagesCellArray = cell(stacksNumber, 1);
    % ROI Selection
    figure('Name', 'Please select the background area')
    M = mat2gray(mean(imagesCellArray{1}(:,:,5:17), 3));
    imshow(M, [min(M(M>0)) .75]);
    r1 = drawrectangle('Label','','Color',[1 0 0]);
    pos = fix(r1.Vertices);
    close
    
    avgs = ones(1, stacksNumber);
    cont = 1;
    for i = 1: stacksNumber
        for j = 1: size(imagesCellArray{i}, 3)
            bkgArea = imagesCellArray{i}(pos(1,2):pos(2,2),pos(1,1):pos(3,1),j);
            bkgMean = mean(bkgArea(:));
            avgs(cont) = bkgMean;
            clearImagesCellArray{i}(:,:,j) = imagesCellArray{i}(:,:,j) - bkgMean;
        
            cont = cont + 1;
        end
    end
    
%     printBackgroundSubstr(avgs, stacksNumber, imagesCellArray, clearImagesCellArray)
end

function printBackgroundSubstr(avgs, stacksNumber, imagesCellArray, clearImagesCellArray)
    figure
    subplot(3,1,1)
    plot(avgs, '.')
    title("Background Value per Frame (In area selected by user)")
    subplot(3,1,2)
    avgs = ones(1, 1800);
    cont = 1;
    for i = 1: stacksNumber
        for j = 1: size(imagesCellArray{i}, 3)
            avgs(cont) = mean(imagesCellArray{i}(:,:,j), "all");
            cont = cont + 1;
        end
    end
    plot(avgs)
    title("Mean Value per Frame of F Images")
    subplot(3,1,3)
    avgs = ones(1, 1800);
    cont = 1;
    for i = 1: stacksNumber
        for j = 1: size(clearImagesCellArray{i}, 3)
            avgs(cont) = mean(clearImagesCellArray{i}(:,:,j), "all");
            cont = cont + 1;
        end
    end
    plot(avgs)
    title("Mean Value per Frame of F Images without Background")
end


