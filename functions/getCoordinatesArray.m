function COORDINATES = getCoordinatesArray(file, urlPath)
%GETCOORDINATESARRAY 
%   GETCOORDINATESARRAY: This function loads the vector of (x,y,frame) coordinates from a
%   txt file
%   COORDINATES = getCoordinatesArray(file, urlPath)
%   returns a N x 3 array, which the x, y and frame coordinate of N ROIs.
%   
%   FILE:    String containing the name of the txt file
%
%   URLPATH: String containing the absolute or relative URL of the
%   txt file location.

    disp("Loading coordinates file... " + file)
    fileId = fopen(urlPath + file);
    tline = fgetl(fileId);
    COORDINATES = zeros(1,3);
    cont = 1;
    while ischar(tline)
        coord = reshape(string(split(tline)), 1, []);
        COORDINATES(cont,:) = coord;
        tline = fgetl(fileId);
        cont = cont + 1;
    end
    fclose(fileId);
end

