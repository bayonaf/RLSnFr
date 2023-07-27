function [selectedSynapses, allCoGs] = getSelectedSynapses(synapseClusters, radius)
%GETSELECTEDSYNAPSES Summary of this function goes here
    
    clusters = unique(synapseClusters(:,1));
    allCoGs = zeros(numel(clusters), 6);
    trials = cell(1);
    for i = 1 : numel(clusters)
        cluster = synapseClusters(synapseClusters(:,1) == clusters(i), :);
        allCoGs(i,1) = clusters(i);
        allCoGs(i,2:5) = [mean(cluster(:,2)), mean(cluster(:,3)), ...
            mean(cluster(:,4)), mean(cluster(:,5))]; 
        trials{i,1} = clusters(i);
        trials{i,2} = unique(cluster(:,6));
    end
    
    % When the distances between synapses is bigger than the radius, is marked as 0
    for i = 1 : size(allCoGs,1)
        flag = false;
        for j = 1 : size(allCoGs,1)
            if i ~= j
                points = [allCoGs(i,2:3) ; allCoGs(j,2:3)];
                dist = pdist(points, 'euclidean');
                if dist <= radius
                    flag = true;
                end
            end
        end
        allCoGs(i,6) = flag;
    end
    
    % The synapses marked with 0 are the isolated ones
    selectedSynapses.CoG = allCoGs(allCoGs(:,6) == 0, :);
    newTrials = cell(1);
    for s = 1 : size(selectedSynapses.CoG,1)
        for t = 1 : size(trials, 1)
            if trials{t,1} == selectedSynapses.CoG(s,1)
                newTrials{s,1} = trials{t,1};
                newTrials{s,2} = trials{t,2};
            end
        end
    end
    selectedSynapses.trials = newTrials;
end

