function eventsPoints = getEvents(varargin)
%GETRELEASECLUSTERS Summary of this function goes here

    switch varargin{1}
        case 0 % triggered
            stacksNumber = length(varargin{2});
            eventsPoints = cell(1, stacksNumber);
    
            for i = 1: stacksNumber
                [row,col,frame] = ind2sub(size(varargin{2}{i}),find(varargin{2}{i} > varargin{3}));
                if ~isempty(row) && ~isempty(col) && ~isempty(frame)
                    tempData = [col, row, frame];
                    tempData(tempData(:, 3) < varargin{4} - 2, :) = []; % Delete all signals before frame 20
                    tempData(tempData(:, 3) > varargin{4} + 2, :) = []; % Delete all signals after frame 22
                    eventsPoints{i}(:,2:4) = tempData;
                    for c = 1 : size(eventsPoints{i}, 1)
                        eventsPoints{i}(c,5) = ...
                            varargin{2}{i}(eventsPoints{i}(c,3), eventsPoints{i}(c,2), eventsPoints{i}(c,4));
                        eventsPoints{i}(c,6) = i;
                    end
                else
                    eventsPoints{i} = [];
                end
            end
        case 1 % spontaneous
            [row,col,frame] = ind2sub(size(varargin{2}),find(varargin{2} > varargin{3}));
            if ~isempty(row) && ~isempty(col) && ~isempty(frame)
                tempData = [col, row, frame];
                eventsPoints(:,2:4) = tempData;
                for c = 1 : size(eventsPoints, 1)
                    eventsPoints(c,5) = ...
                        varargin{2}(eventsPoints(c,3), eventsPoints(c,2), eventsPoints(c,4));
                end
                eventsPoints(:,6) = NaN;
            else
                eventsPoints = [];
            end
            eventsPoints = {eventsPoints};
    end
end

