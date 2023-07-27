function [clusters, releaseClusters] = getClusterIndexes(releaseClusters, epsilon, minPts, frameFactor, clusteringCase)
    
    numStacks = size(releaseClusters, 2);
    switch clusteringCase
        case 'Release Clusters'
            clusterCont = 0;
            for i = 1: numStacks
                if ~isempty(releaseClusters{i})
                    % Fram multipliyed by factor to make epsilon in time
                    % bigger than spatial epsilon
                    releaseClusters{i}(:,4) = releaseClusters{i}(:,4) * frameFactor; 
                    % DBScan clustering
                    idxDBSCAN = dbscan(releaseClusters{i}(:,2:4), epsilon, minPts);
                    % Frame default value is restored
                    releaseClusters{i}(:,4) = releaseClusters{i}(:,4) / frameFactor; 
                    releaseClusters{i}(:,1) = idxDBSCAN;
                    releaseClusters{i}(:,1) = releaseClusters{i}(:,1) + clusterCont;
                    clusterCont = max(releaseClusters{i}(:,1));
                end
            end
            clusters = releaseClusters;
        case 'Synapse Clusters'
            CoGs = getCentersOfGravity(releaseClusters);
            % All centers of gravity concatenated in a single matrix (across trials)
            clusters = cat(1, CoGs{:}); 
            % DBScan execution over x,y of CoGs
            clusters(:,1) = dbscan(clusters(:,2:3), epsilon, minPts); 
            releaseClusters = propagateIndexes(clusters(:,1), CoGs, releaseClusters);
    end
end

% Propagate the new clusters number to the release clusters after
% clustering synapse clusters
function releaseClusters = propagateIndexes(idx, CoGs, releaseClusters)
    cuttingCont = 1;
    allIdx = cell(0);
    for i = 1 : length(CoGs)
        s = size(CoGs{i},1);
        allIdx{i} = [CoGs{i}(:,1), idx(cuttingCont : cuttingCont + s -1, :)];
        cuttingCont = cuttingCont + s;
    end
    for i = 1 : length(releaseClusters)
        for rc = 1 : size(releaseClusters{i}, 1)
            ind = allIdx{i}(:,1) == releaseClusters{i}(rc,1);
            releaseClusters{i}(rc,1) = allIdx{i}(ind,2);
        end
    end
end
