function exportResultsToCSV(releaseClusters, selectedSynapses, ROIs, unbleachedROIs, unbleachedNormalizedROIs, ROIsFitCoeffs, selpath, numTrials)
% EXPORTRESULTSTOCSV(): This function exports a the CVS and XLS files
% containing the selected synapses traces and information
%   
%   EXPORTRESULTSTOCSV(selectedSynapses, unbleachedROIs, ROIsFitCoeffs, originalROIs, selpath)
%   
%   SELECTEDSYNAPSES:    Cell array containing all the selected synapses  
%   UNBLEACHEDROIS:      Cell array containing all the synapses DF/F
%   unbleached traces
%   ROIsFitCoeffs:       Cell array containing the fitting coefficients for
%   each trace per trial after unbleaching and normalization
%   ORIGINALROIS:        Cell array containing all the synapse traces per
%   trial before normalization
%   SELPATH:             Selected URL to save the file
    
    if ~exist(selpath  + "\3. Normalized Synapses\", 'dir')
        mkdir(selpath + "\3. Normalized Synapses\");
        mkdir(selpath + "\2. Corrected Synapses\");
        mkdir(selpath + "\1. Raw Synapses\");
    end
    synapseNrs = selectedSynapses.CoG(:,1);

    % ===== Synapses Correction and Normalization Fitting Coefficients
    titles = ["Correction Type", "Baselien Decay %", "A", "B", "tau_decay", "tau_bleach","y_offset", "M", "pseudo_B"];
    leg    = ["A: Initial Absolute Amplitude"; "B: Relative Peak Amplitude"; "tau_decay: Signal Decay Rate"; ...
        "tau_bleach: Baseline Exp. Decay Rate"; "yoffset: Bleaching coefficient"; "tau_bleach: Baseline Linear Decay Rate";...
        "pseudo_B: Baseline Relative pseud Amplitude"];
    Fs = ["Correction Type 1: ((1 + (B * exp (-(x-xresp) ./ tau_decay))) .* (x>=xresp) + (x<xresp)) * A"; ...
        "Correction Type 2: ((1 + (B * exp (-(x-xresp) ./ tau_decay))) .* (x>=xresp) + (x<xresp)) * A .* ((M * (x-(xbegin-1))) + 1)"; ...
        "Correction Type 3: ((1 + (B * exp (-(x-xresp) ./ tau_decay))) .* (x>=xresp) + (x<xresp)) * A .* ((exp (-(x-(xbegin-1)) ./ tau_bleach)) .* (1-yoffset) + yoffset)";...
        "";"";"";""];
    leg = [leg, strings(7,1), Fs];
    filename = selpath + "Correction-Normalization Coefficients.xlsx";
    for R = 1 : size(ROIsFitCoeffs, 1)
        writematrix(leg, filename, "Sheet", "Coeffs " + string(synapseNrs(R)), "Range", "A1");
        writematrix("                           ", filename, "Sheet", "Coeffs " + string(synapseNrs(R)), "Range", "A8")
        writematrix(titles, filename, "Sheet", "Coeffs " + string(synapseNrs(R)), "Range", "A9");
        writematrix(round(ROIsFitCoeffs{R},3), filename, "Sheet", "Coeffs " + string(synapseNrs(R)), "Range", "A10");
        writematrix("                           ", filename, "Sheet", "Coeffs " + string(synapseNrs(R)), 'WriteMode', 'append')
    end


%     % ===== Release Clusters of Selected Synapses
%     titles = ["Synapse Nr.", "X", "Y", "Frame", "Dec. Peak Amp.", "Trial"];
%     filename = selpath + "Release Clusters.xlsx";
%     clusters = cat(1, releaseClusters{:}); 
%     strs = strings(1,length(titles));
%     for i = 1 : length(strs); strs(i) = "                   ";end
%     for i = 1 : length(synapseNrs)
%         cluster = round(clusters(clusters(:, 1) == synapseNrs(i),:), 1);
%         writematrix(titles, filename, "Sheet", "Synapse " + string(synapseNrs(i)))
%         writematrix(cluster, filename, "Sheet", "Synapse " + string(synapseNrs(i)), 'WriteMode', 'append')
%         writematrix(strs, filename, "Sheet", "Synapse " + string(synapseNrs(i)), 'WriteMode', 'append')
%     end

    % ===== Synapses Centers of Gravity, Areas, and Release Trials
    paramsMap = getParametersMap(); 
    filename = selpath + "Synapse Positions.xlsx";
    titles = {"SYNAPSES CoG"; 
        ["Synapse Nr.", "CoG X pos.", "CoG Y pos.",...
        "CoG Time", "Syn. Max. Area", "Syn. Min. Area"];
        ["", "px", "px", "frame", "nm^2", "nm^2"]};
    trialsTitles = strings(1, numTrials + 1);
    for i = 1 : numTrials 
        trialsTitles(i) = "Trial " + string(i);
    end
    titles{2} = [titles{2}, "", trialsTitles];
    strs = strings(1,length(titles{2}));
    for i = 1 : length(strs); strs(i) = "                   "; end
    synapseInfo = [round(selectedSynapses.CoG(:,1:4),2), round(selectedSynapses.CoG(:,6:7),4)];
    trialsMatrix = NaN(length(selectedSynapses.trials), numTrials);
    for i = 1 : length(selectedSynapses.trials)
        trialsMatrix(i, selectedSynapses.trials{i,2}) = 1;
    end
    data = [synapseInfo, NaN(size(synapseInfo,1),1), trialsMatrix];
    data = num2cell(data);
    temp = cellfun(@ismissing,data);
    data(temp) = {[]};
    writematrix([titles{1}], filename, "Sheet", "CoG - Areas", "Range", "A1")
    writematrix(sprintf("Px = %.2fum)", paramsMap('pixelFactorNm')/1000), filename, "Sheet", "CoG - Areas", "Range", "A2")
    writematrix(titles{2}, filename, "Sheet", "CoG - Areas", "Range", "A3")
    writematrix(titles{3}, filename, "Sheet", "CoG - Areas", "Range", "A4")
    writecell(data, filename, "Sheet", "CoG - Areas", "Range", "A5")
    writematrix(strs, filename, "Sheet", "CoG - Areas", 'WriteMode', 'append')

    
    % ===== Normalized Traces
    % For users reading
    titles = {"Synapse Nr."; "Trace Type"};
    for i = 3 : size(unbleachedNormalizedROIs{1}, 2) + 1
        titles{i,1} = "Frame " + string(i-2);
    end
    titles{end+1,1} = "Initial Abs. Amplitude";
    titles{end+1,1} = "Rel. Response Amplitude"; 
    titles{end+1,1} = "Decay Rate";
    s = size(unbleachedNormalizedROIs{1}(3:end,:)');
    synapseNr = NaN(1,s(2) + 5);
    data = NaN(s(1)+4 ,1);
    for R = 1 : size(unbleachedNormalizedROIs, 1)
        synapseNr(1) = synapseNrs(R); 
        D = round(unbleachedNormalizedROIs{R}(3:end,:)',4);
        D = cat(2, D, NaN(s(1) ,1));
        D = cat(2, D, round(unbleachedNormalizedROIs{R}(1:2,:)',4));
        D = cat(2, D, NaN(s(1) ,2));
        C = cat(2, round(ROIsFitCoeffs{R}(:,3:5)',2), NaN(3,5));
        D = cat(1, synapseNr, D, C );
        data = cat(2, data, D);
    end
    data = num2cell(data);
    temp = cellfun(@ismissing,data);
    data(temp) = {[]};
    data(:,1) = titles;
    filename = selpath + "All Normalized Synapses.xlsx";
    writecell(data, filename)
    writematrix(" ", filename,'WriteMode','append')
    Leg = ["Legend:"; "Type 1: Response trace"; "Type 2: Failure trace"];
    writematrix(Leg, filename,'WriteMode','append')

    % For processing
    s = size(unbleachedNormalizedROIs{1}(3:end,:)');
    data = [];
    for R = 1 : size(unbleachedNormalizedROIs, 1)
        titles = strings(1, numTrials);
        for i = 1 : numTrials
            titles(i) = sprintf("Norm%0i-T%0i",[synapseNrs(R), i]);
        end
        D = round(unbleachedNormalizedROIs{R}(3:end,2:end)',4);
        data = cat(1, data, titles, D, NaN(1, s(2)));
    end
    data = num2cell(data);
    temp = cellfun(@ismissing,data);
    data(temp) = {[]};
    filename = selpath + "All Normalized Synapses.csv";
    writecell(data, filename)


    % ===== Corrected Traces
    % For users reading
    titles = {"Synapse Nr."; "Trace Type"};
    for i = 3 : size(unbleachedROIs{1}, 2) + 1
        titles{i,1} = "Frame " + string(i-2);
    end
    titles{end+1,1} = "Initial Abs. Amplitude";
    titles{end+1,1} = "Rel. Response Amplitude"; 
    titles{end+1,1} = "Decay Rate";
    s = size(unbleachedROIs{1}(3:end,:)');
    synapseNr = NaN(1,s(2) + 5);
    data = NaN(s(1)+4 ,1);
    for R = 1 : size(unbleachedROIs, 1)
        synapseNr(1) = synapseNrs(R); 
        D = round(unbleachedROIs{R}(3:end,:)',1);
        D = cat(2, D, NaN(s(1) ,1));
        D = cat(2, D, round(unbleachedROIs{R}(1:2,:)',1));
        D = cat(2, D, NaN(s(1) ,2));
        C = cat(2, round(ROIsFitCoeffs{R}(:,3:5)',2), NaN(3,5));
        D = cat(1, synapseNr, D, C );
        data = cat(2, data, D);
    end
    data = num2cell(data);
    temp = cellfun(@ismissing,data);
    data(temp) = {[]};
    data(:,1) = titles;
    filename = selpath + "All Corrected Synapses.xlsx";
    writecell(data, filename)
    writematrix(" ", filename,'WriteMode','append')
    Leg = ["Legend:"; "Type 1: Response trace"; "Type 2: Failure trace"];
    writematrix(Leg, filename,'WriteMode','append')

    % For processing
    s = size(unbleachedROIs{1}(3:end,:)');
    data = [];
    for R = 1 : size(unbleachedROIs, 1)
        titles = strings(1, numTrials);
        for i = 1 : numTrials
            titles(i) = sprintf("Corr%0i-T%0i",[synapseNrs(R), i]);
        end
        D = round(unbleachedROIs{R}(3:end,2:end)',1);
        data = cat(1, data, titles, D, NaN(1, s(2)));
    end
    data = num2cell(data);
    temp = cellfun(@ismissing,data);
    data(temp) = {[]};
    filename = selpath + "All Corrected Synapses.csv";
    writecell(data, filename)


    % ===== Raw Traces
    % For users reading
    titles = {"Synapse Nr."; "Trace Type"};
    for i = 3 : size(ROIs{1}, 2) + 1
        titles{i,1} = "Frame " + string(i-2);
    end
    titles{end+1,1} = "Initial Abs. Amplitude";
    titles{end+1,1} = "Rel. Response Amplitude"; 
    titles{end+1,1} = "Decay Rate";
    s = size(ROIs{1}');
    synapseNr = NaN(1,s(2) + 2);
    data = NaN(s(1)+4 ,1);
    for R = 1 : size(ROIs, 1)
        synapseNr(1) = synapseNrs(R); 
        D = round(single(ROIs{R})',1);
        D = cat(2, D, NaN(s(1) ,2));
        C = cat(2, round(ROIsFitCoeffs{R}(:,3:5)',2), NaN(3,2));
        D = cat(1, synapseNr, D, C );
        data = cat(2, data, D);
    end
    data = num2cell(data);
    temp = cellfun(@ismissing,data);
    data(temp) = {[]};
    data(:,1) = titles;
    filename = selpath + "All Raw Synapses.xlsx";
    writecell(data, filename)
    writematrix(" ", filename,'WriteMode','append')
    Leg = ["Legend:"; "Type 1: Response trace"; "Type 2: Failure trace"];
    writematrix(Leg, filename,'WriteMode','append')

    % For processing
    s = size(ROIs{1}');
    data = [];
    for R = 1 : size(ROIs, 1)
        titles = strings(1, numTrials);
        for i = 1 : numTrials
            titles(i) = sprintf("Raw%0i-T%0i",[synapseNrs(R), i]);
        end
        D = round(single(ROIs{R}(:,2:end))',1);
        data = cat(1, data, titles, D, NaN(1, s(2)));
    end
    data = num2cell(data);
    temp = cellfun(@ismissing,data);
    data(temp) = {[]};
    filename = selpath + "All Raw Synapses.csv";
    writecell(data, filename)


    % ===== Export ROI per file
    trialsTitles = strings(1, size(unbleachedROIs{1}, 1));
    for i = 1 : size(unbleachedROIs{1}, 1) - 2
        trialsTitles(i) = "Trial " + string(i);
    end
    trialsTitles(end-1) = "All Signals Avg";
    trialsTitles(end) = "Resp Signals Avg";

    parfor R = 1 : size(unbleachedNormalizedROIs, 1)
        exportROItoCSV(unbleachedNormalizedROIs{R}', selpath, 'Normalized', trialsTitles, synapseNrs(R))
    end
    parfor R = 1 : size(unbleachedROIs, 1)
        exportROItoCSV(unbleachedROIs{R}', selpath, 'Corrected', trialsTitles, synapseNrs(R))
    end
    parfor R = 1 : size(unbleachedROIs, 1)
        exportROItoCSV(single(ROIs{R}'), selpath, 'Raw', trialsTitles, synapseNrs(R))
    end
end
    

function exportROItoCSV(ROI, selpath, roiType, trialsTitles, synapseNr)
    titles = ["SYNAPSES TRACES"; "Synapse nr."];
    
    switch roiType
        case 'Normalized'
            selpath = selpath + "\3. Normalized Synapses\";
            filename = selpath + "NormROI" + synapseNr + ".csv";
            ROI = round(ROI, 4);
            temp = ROI(:,1:2);
            ROI(:,1:end-2) = ROI(:,3:end);
            ROI(:,end-1:end) = temp;
        case 'Corrected'
            selpath = selpath + "\2. Corrected Synapses\";
            filename = selpath + "CorrROI" + synapseNr + ".csv";
            ROI = round(ROI, 1);
            temp = ROI(:,1:2);
            ROI(:,1:end-2) = ROI(:,3:end);
            ROI(:,end-1:end) = temp;
        case 'Raw'
            selpath = selpath + "\1. Raw Synapses\";
            filename = selpath + "RawROI" + synapseNr + ".csv";
            trialsTitles(end-1:end) = [];
            ROI = round(ROI, 1);
    end
    
    writematrix(titles(1), filename)
    writematrix("", filename,'WriteMode','append')
    writematrix(titles(2), filename,'WriteMode','append')
    writematrix(synapseNr, filename,'WriteMode','append')
    writematrix("", filename,'WriteMode','append')
    writematrix(trialsTitles, filename,'WriteMode','append')
    writematrix(ROI, filename,'WriteMode','append')
end    