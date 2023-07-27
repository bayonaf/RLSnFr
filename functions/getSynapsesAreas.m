function [selectedSynapses, nrTrialsTriggered] = getSynapsesAreas(synapseClusters, selectedSynapses, pixelFactorNm)
%GETSYNAPSESAREAS Summary of this function goes here
    clusterNr = selectedSynapses.CoG(:,1);
    pixelArea = (pixelFactorNm/1000) ^ 2;
    nrTrialsTriggered = zeros(length(clusterNr), 1);
    for i = 1 : length(clusterNr)
        clus = synapseClusters(synapseClusters(:,1) == clusterNr(i), :);
        nrTrialsTriggered(i) = numel(unique(clus(:,6)));
        try
            [~,av] = convhull(clus(:, 2:3));
            selectedSynapses.CoG(i,6) = av * pixelArea;
        catch
            dist = pdist(clus(:,2:3));
            if dist == 1
                selectedSynapses.CoG(i,6) = pixelArea * 0.2;
            else 
                selectedSynapses.CoG(i,6) = pixelArea * 0.1;
            end
        end
    end
end

