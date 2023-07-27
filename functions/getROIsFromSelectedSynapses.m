function ROIs = getROIsFromSelectedSynapses(filteredStacks, selectedSynapses)
%GETROISFROMSELECTEDSYNAPSES Summary of this function goes here
    numStacks = length(filteredStacks);
    ROIs = cell(size(selectedSynapses.CoG, 1),1);
    for c = 1 : size(selectedSynapses.CoG, 1)
        trials = unique(selectedSynapses.trials{c,2});
        allProfiles = [];
        selectedSynapses.CoG(c,2:3) = round(selectedSynapses.CoG(c,2:3));
        for i = 1 : numStacks
            timeProfile = reshape(filteredStacks{i}(selectedSynapses.CoG(c,3),...
                selectedSynapses.CoG(c,2), :),1,[]);
            if ismember(i, trials)
                allProfiles = cat(1, allProfiles, [1, timeProfile]);
            else
                allProfiles = cat(1, allProfiles, [0, timeProfile]);
            end
        end
        ROIs{c} = allProfiles;
    end
end

