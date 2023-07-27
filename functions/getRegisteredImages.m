function [regImagesCellArray, regCoeffs] = getRegisteredImages(imagesCellArray, stacksNumber, RegistrationMode)
%GETREGISTEREDIMAGES Summary of this function goes here
    
%     tic
%     disp("Registering images... ")
    [optimizer,metric] = imregconfig(RegistrationMode);

    fixed = mean(imagesCellArray{1}(:,:,4:16),3);
    [regImagesCellArray, regInfo] = deal(cell(1));
    parfor j = 1 : stacksNumber
        registered = zeros(size(imagesCellArray{j}), 'single');
        reg = cell(1);
        for i = 1 : size(imagesCellArray{j}, 3)
            moving = imagesCellArray{j}(:,:,i);
            tform = imregtform(moving,fixed,'translation',optimizer,metric);
            registered(:,:,i) = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
            reg{i} = tform;
        end
        regImagesCellArray{j} = registered;
        regInfo{j} = reg;
    end
    regCoeffs = zeros(2, stacksNumber*size(imagesCellArray{1},3));
    cont = 1;
    for i = 1 : stacksNumber
        for j = 1 : size(imagesCellArray{1},3)
            regCoeffs(1, cont) = regInfo{i}{j}.T(3,1);
            regCoeffs(2, cont) = regInfo{i}{j}.T(3,2);
            cont = cont + 1;
        end
    end
%     toc
end

