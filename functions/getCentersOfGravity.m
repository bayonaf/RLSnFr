function CoGs = getCentersOfGravity(releaseClusters)
%GETCENTERSOFGRAVITY Summary of this function goes here
    
    numStacks = size(releaseClusters, 2);
    CoGs = cell(1, numStacks);
    
    for i = 1 : numStacks
        CoGs{i} = []; % Empty vector for all clustered CoG
        if ~isempty(releaseClusters{i})
            uniqueClusters = unique(releaseClusters{i}(:, 1)); % Unique cluster nr. in trial i
            clusters = zeros(numel(uniqueClusters), 5); % Memory allocation
            for j = 1 : numel(uniqueClusters)
                % Decon points info for cluster nr. j in trial i
                elems = releaseClusters{i}(releaseClusters{i}(:, 1) == uniqueClusters(j),:);
                % Center of gravity by averages
                clusters(j,2:6) = [mean(elems(:,2)), mean(elems(:,3)), ...
                    mean(elems(:,4)), mean(elems(:,5)), elems(1,6)]; 
                clusters(j,1) = uniqueClusters(j); % Cluster nr.
            end
            CoGs{i} = cat(1, CoGs{i}, clusters); % All CoG per trial saved
        end
    end
end

