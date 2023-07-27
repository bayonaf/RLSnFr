function [selectedSynapses, deletedSynapses] = getReliableSynapses(selectedSynapses, clusterAreas, synapseAreas, clustersLabels, synapseLabels)
%GETPROPPERSYNAPSES Summary of this function goes here
%   Detailed explanation goes here
    
    selectedSynapses.CoG(:,6) = clusterAreas(selectedSynapses.CoG(:,1));
    selectedSynapses.CoG(:,7) = synapseAreas(selectedSynapses.CoG(:,1));
    keptSynapses = selectedSynapses;
    releableLabels = logical(clustersLabels) | logical(synapseLabels);
    idx = ismember(keptSynapses.CoG(:,1), find(releableLabels == 1));
    keptSynapses.CoG(idx,:) = [];
    keptSynapses.trials(idx,:) = [];

    % Temporal
    idx = ismember(selectedSynapses.CoG(:,1), find(releableLabels == 0));
    deletedSynapses = selectedSynapses;
    deletedSynapses.CoG(idx,:) = [];
    deletedSynapses.trials(idx,:) = [];
    % Temporal
end

