function M = getParametersMap()
%GETPARAMETERSMAP Summary of this function goes here
    if isunix || ismac
        fileID = fopen('data/key_value.m','r');
    elseif ispc
        fileID = fopen('data\key_value.m','r');
    end
    allLines = string(1);
    pos = 1;
    while true
        readLine = fgetl(fileID);
        if readLine == -1
            break
        end
        if(~strcmp(readLine(1),'%'))
            allLines(pos,1) = readLine;
            pos = pos + 1;
        end
    end
    fclose(fileID);

    splited = split(allLines,'=');
    M = containers.Map(splited(:,1), str2double(splited(:,2)));
end

