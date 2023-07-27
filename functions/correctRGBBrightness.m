function imagesToVisualize = correctRGBBrightness(imagesToVisualize, stacksNumber)
    % Getting Maximum to Normalize
    M = max(max(imagesToVisualize.Original, imagesToVisualize.Registered), [], 'all');

    % Normalizing Original
    imagesToVisualize.Original = imagesToVisualize.Original / M;
    % Getting Minimum of Original
    [~, ind] = min(imagesToVisualize.Original,[],"all",'linear');
    [row, col, ~] = ind2sub(size(imagesToVisualize.Original), ind);
    try
        m = imagesToVisualize.Original(row-3:row+3, col-3:col+3, :);
        m = [mean(m(:,:,1),"all"), mean(m(:,:,2),'all'), mean(m(:,:,3),"all")];
    catch
        m = imagesToVisualize.Original(row, col,:);
        m = [m(1), m(2), m(3)];
    end
    % Substracting Original Offset
    imagesToVisualize.Original(:,:,1) = imagesToVisualize.Original(:,:,1) - m(1);
    imagesToVisualize.Original(:,:,2) = imagesToVisualize.Original(:,:,2) - m(2);
    imagesToVisualize.Original(:,:,3) = imagesToVisualize.Original(:,:,3) - m(3);
    
    % Normalizing Registered
    imagesToVisualize.Registered = imagesToVisualize.Registered / M;
    % Substracting Registered Offset
    imagesToVisualize.Registered(:,:,1) = imagesToVisualize.Registered(:,:,1) - m(1);
    imagesToVisualize.Registered(:,:,2) = imagesToVisualize.Registered(:,:,2) - m(2);
    imagesToVisualize.Registered(:,:,3) = imagesToVisualize.Registered(:,:,3) - m(3);

    % Normalizing Template
    imagesToVisualize.Template = (imagesToVisualize.Template / M) / stacksNumber;
    imagesToVisualize.Template = imagesToVisualize.Template  - mean(m);
end

