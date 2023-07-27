function IMAGESMATRIXCELL = getMatrixFromImages(imageFile, imageUrlPath, isParallel)
% GETMATRIXFROMIMAGES: This function reads the tiff images stack from a url
%   that comes like function argument, and returns a matlab 3D matrix with
%   images information.
%   IMAGESMATRIXCELL = getMatrixFromImages(imagesFile, imagesUrlPath, isParallel) 
%   returns a num_images x 1 cell array, which each element is a 3 dimensional matrix
%   with the numerical values of the original TIFF stacks
%   
%   IMAGEFILE:    String containing the name of the TIFF file
%
%   IMAGEURLPATH: String containing the absolute or relative URL of the
%   stack TIFF file location.
%
%   ISPARALLEL:   Boolean that defines if the images loading must be executed 
%   in parallel or 1 CPU  

    if isParallel
        info = imfinfo(imageUrlPath + imageFile);
        numImages = numel(info); % Number of images in the tiff file
        IMAGESMATRIXCELL = cell(numImages, 1); % Cell matrix of dimension num_images x 1
        
        parfor (i = 1 : numImages) % Parallel for loop
            A = imread(imageUrlPath + imageFile, i);
            IMAGESMATRIXCELL{i} = A; % Insert the images in the matrix
        end
        
    else
        info = imfinfo(imageUrlPath + imageFile);
        numImages = numel(info); % Number of images in the tiff file
        IMAGESMATRIXCELL = cell(numImages, 1); % Cell matrix of dimension num_images x 1
        
        for i = 1 : numImages
            A = imread(imageUrlPath + imageFile, i);
            IMAGESMATRIXCELL{i} = A; % Insert the images in the matrix
        end
    end
    
end


