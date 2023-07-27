% The purpose of this test is to measure the time and memory usage of the 
% tool during each stage of the analysis. Will be assumed a standard PC
% architecture with 4 cores, and will be analysed an experiment consisting
% of 4 trials of size 512 x 512 x 30 to evaluate the machine requirements. 
% The goal is to assess the tool's performance and resource needs 

% Loading Parameters
paramsMap = getParametersMap();
RLIter = paramsMap('RLIterations');
triggerTime = paramsMap('triggerTime');
baselineEndFrame = paramsMap('baselineEndFrame');
initFrames = paramsMap('initFrames');
regMode = paramsMap('regMode'); if regMode == 0; regMode = 'monomodal'; else; regMode = 'multimodal'; end

% General Counter
gCont = 1;
% Statistics
timeConsumption = zeros(0);
memConsumption  = zeros(0);

% ===== Load Experiment Data
t = tic;
disp("1. Loading Images...")
stacksNumber = 4;
fileNames = ["0.tif", "1.tif", "2.tif", "3.tif"];
imagesCellArray = cell(stacksNumber, 1); % Cell array containing all the 3D stacks loaded
for i = 1: stacksNumber
    imageMatrixCell = getMatrixFromImages(fileNames(i), "data/TestFiles/", false);
    imagesCellArray{i, 1} = cat(3,imageMatrixCell{:});
end
imagesToVisualize.Original = getShiftedCompositeImage(imagesCellArray);
imagesToVisualize.Template = mean(imagesCellArray{1}(:,:,initFrames-1:triggerTime), 3);
M = memory; 
memConsumption(gCont) = M.MemUsedMATLAB;
timeConsumption(gCont)  = toc(t);
gCont = gCont + 1;
toc(t)
clear fileNames imageMatrixCell i M


%% ===== Images Registration
t = tic;
disp("2. Registering images... ")
[regImagesCellArray, regCoeffs] = getRegisteredImages(imagesCellArray, stacksNumber, "monomodal");
imagesToVisualize.Registered = getShiftedCompositeImage(regImagesCellArray);
imagesToVisualize = correctRGBBrightness(imagesToVisualize, stacksNumber);
M = memory; 
memConsumption(gCont) = M.MemUsedMATLAB;
timeConsumption(gCont)  = toc(t);
gCont = gCont + 1;
toc(t)
clear M t

%% ===== Bleaching Correction
t = tic;
disp("3. Unbleaching Images...")
[unbleachedImagesCellArray, UnbleachingResults.meanPercentDiff, ...
    UnbleachingResults.avgBefore, UnbleachingResults.avgAfter, UnbleachingResults.avgWhole,...
    UnbleachingResults.fits, UnbleachingResults.allCoeffs]...
    = getImagesWithoutBleaching(regImagesCellArray, initFrames-1, triggerTime);
M = memory; 
memConsumption(gCont) = M.MemUsedMATLAB;
timeConsumption(gCont)  = toc(t);
gCont = gCont + 1;
toc(t)
clear M t


%% ===== DF Images
t = tic;
disp("4. Getting ΔF Images...")
[deltaFStacks, deltaFOverFstacks, baselinesRegCellArray] = deal(cell(stacksNumber,1));
for i = 1: stacksNumber
    [deltaFStacks{i, 1}, deltaFOverFstacks{i, 1}, baselinesRegCellArray{i, 1}] =...
        getDFMatrix(0, unbleachedImagesCellArray{i}, initFrames-1, baselineEndFrame);
end
M = memory; 
memConsumption(gCont) = M.MemUsedMATLAB;
timeConsumption(gCont)  = toc(t);
gCont = gCont + 1;
toc(t)
clear M t i


%% ===== Deconvolution
t = tic;
disp("5. Deconvolution Started... ")
load("data/Triggered Kernel.mat", "defaultKernel")
deconvStacks = cell(stacksNumber, 1);
parfor i = 1: stacksNumber
    deconvStacks{i, 1} = deconvlucy(deltaFStacks{i}, defaultKernel, RLIter);
end
M = memory; 
memConsumption(gCont) = M.MemUsedMATLAB;
timeConsumption(gCont)  = toc(t);
gCont = gCont + 1;
toc(t)
clear M t i


%% ===== RESULTS ANALYSIS
experimentType = paramsMap('experimentType');
threshFactor = paramsMap('threshFactor');
epsilonRCnm = paramsMap('epsilonRCnm');
epsilonSCnm = paramsMap('epsilonSCnm');
gaussSigmaSAnm = paramsMap('gaussSigmaSA');
radiusSCnm = paramsMap('radiusSCnm');
pixelFactorNm = paramsMap('pixelFactorNm');
minPts = paramsMap('minPts');
epsilonRC = epsilonRCnm/pixelFactorNm;
epsilonSC = epsilonSCnm/pixelFactorNm;
gaussSigmaSA = gaussSigmaSAnm/pixelFactorNm;
radiusSC = radiusSCnm/pixelFactorNm;
pseudoTriggerTime = paramsMap("pseudoTriggerTime");

%% ===== Get Threshold, potential release events, and Release Clusters
t = tic;
disp("6. Getting Potential Release Events... ")
filteredStacks = cell(1); 
parfor i = 1 : stacksNumber
    filteredStacks{i} = imgaussfilt(deltaFStacks{i}, gaussSigmaSA);
end
allPos = zeros(1, stacksNumber);
for i = 1 : stacksNumber
    temp = deconvStacks{i}(:,:,initFrames:baselineEndFrame);
    temp = temp(temp > 0);
    allPos(i) = prctile(temp, 95, 'all');
end
threshold = mean(allPos,'omitnan');
ThresholdVal = threshold * threshFactor;
eventsPoints = getEvents(experimentType, deconvStacks, ThresholdVal, triggerTime);
eventsPoints = getEventsWithoutOutliers(experimentType, eventsPoints, regCoeffs, size(deconvStacks{1},1), [], []);
trialsIndices = getCorrectPotentialEvents(filteredStacks, eventsPoints, initFrames, triggerTime, pseudoTriggerTime);
for i = 1: stacksNumber
    eventsPoints{i} = eventsPoints{i}(logical(trialsIndices{i}), :);
end
releaseClusters = getClusterIndexes(eventsPoints, epsilonRC, minPts, 1, 'Release Clusters');
centersOfGravityRC = getCentersOfGravity(releaseClusters);
M = memory; 
memConsumption(gCont) = M.MemUsedMATLAB;
timeConsumption(gCont)  = toc(t);
gCont = gCont + 1;
toc(t)
clear M t i temp allPos trialsIndices

% ===== Get Synapse Clusters
t = tic;
disp("7. Getting Synapses... ")
[synapseClusters, releaseClusters] = getClusterIndexes(releaseClusters, epsilonSC, minPts, [], 'Synapse Clusters');
[clusterAreas, clustersLabels] = getClustersAreas(releaseClusters, pixelFactorNm, pixelFactorNm*3, 'Release Clusters');
[synapseAreas, synapseLabels] = getClustersAreas(synapseClusters, pixelFactorNm, pixelFactorNm*2, 'Synapse Clusters');
[selectedSynapses, SynapsesCoG] = getSelectedSynapses(synapseClusters, radiusSC);
selectedSynapses = getReliableSynapses(selectedSynapses, clusterAreas, synapseAreas, clustersLabels, synapseLabels);
M = memory; 
memConsumption(gCont) = M.MemUsedMATLAB;
timeConsumption(gCont)  = toc(t);
gCont = gCont + 1;
toc(t)
clear M t i clustersLabels


%% ===== Results export
t = tic;
disp("8. Exporting Results... ")
filteredStacks = cell(1); 
[FImages, Backgrs, vertices] = getImagesWithoutBackground(regImagesCellArray);
parfor i = 1 : stacksNumber
    filteredStacks{i} = imgaussfilt(FImages{i}, gaussSigmaSA);
end
ROIs = getROIsFromSelectedSynapses(filteredStacks, selectedSynapses);
[unbleachedROIs, unbleachedNormalizedROIs, ROIsFitCoeffs] = getUnbleachedNormalizedROis(ROIs, initFrames-1, triggerTime);
[unbleachedNormalizedROIs, unbleachedROIs] = getROIsAverages(unbleachedNormalizedROIs, unbleachedROIs);
exportResultsToCSV(releaseClusters, selectedSynapses, ROIs, unbleachedROIs, unbleachedNormalizedROIs, ROIsFitCoeffs, "data\TestFiles\", stacksNumber)
exportMetadata("data\TestFiles\", "PerformanceTest", threshFactor, ThresholdVal, vertices, Backgrs, initFrames-1, paramsMap)
exportImages("data\TestFiles\", regImagesCellArray, deltaFStacks, imagesToVisualize, regCoeffs,...
    SynapsesCoG, selectedSynapses, unbleachedNormalizedROIs, vertices)
exportToTIF(deconvStacks, "data\TestFiles\", "0. Images\Deconv Stacks", "\Dec ")
exportToTIF(deltaFStacks, "data\TestFiles\", "0. Images\DF Stacks", "\DF ")
exportToTIF(deltaFOverFstacks, "data\TestFiles\", "0. Images\DFF Stacks", "\DFF ")
exportToTIF(regImagesCellArray, "data\TestFiles\", "0. Images\Registered Stacks", "\Reg ")
M = memory; 
memConsumption(gCont) = M.MemUsedMATLAB;
timeConsumption(gCont)  = toc(t);
gCont = gCont + 1;
toc(t)
clear M t i

% ===== Exporting Performance Info
exportPerformanceData(memConsumption, timeConsumption);


%% ===== AUX FUNCTIONS =====
function [unbleachedNormalizedROIs, unbleachedROIs] = getROIsAverages(unbleachedNormalizedROIs, unbleachedROIs)
    for i = 1 : size(unbleachedROIs, 1)
        M = mean(unbleachedROIs{i});
        M(1) = 2;
        realSignals = unbleachedROIs{i}(unbleachedROIs{i}(:,1) == 1, :);
        Ms = mean(realSignals, 1, "omitnan");
        Ms(1) = 3;
        unbleachedROIs{i} = cat(1, M, Ms, unbleachedROIs{i});

        M = mean(unbleachedNormalizedROIs{i});
        M(1) = 2;
        realSignals = unbleachedNormalizedROIs{i}(unbleachedNormalizedROIs{i}(:,1) == 1, :);
        Ms = mean(realSignals, 1, 'omitnan');
        Ms(1) = 3;
        unbleachedNormalizedROIs{i} = cat(1, M, Ms, unbleachedNormalizedROIs{i});
    end
end

function exportImages(selpath, registeredStacks, deltaFStacks, imagesToVisualize, regCoeffs,...
    SynapsesCoG, selectedSynapses, unbleachedNormalizedROIs, vertices)
    trilaColors = {[0,0.447,0.741, 0.3]
            [0.85,0.325,0.098,1]};
    % Plot Before and After Registration
    stacksNumber = length(registeredStacks);
    Ticks = linspace(0,stacksNumber+1, stacksNumber+1);
    TicksLabels = string(linspace(0,stacksNumber, stacksNumber+1));
    f = figure('units','normalized','outerposition',[0 0 1 1], 'Visible','off');
    sgtitle("Images Before and After Registration")
    ax(1) = subplot(1,2,1);
    imshow(imagesToVisualize.Original)
    cmap = hsv(stacksNumber);
    colormap(ax(1), cmap)
    c = colorbar;
    c.Ticks = Ticks/max(Ticks);
    c.TickLabels = TicksLabels;
    ylabel(c,"Trials");
    title("Raw Data Composite Image")
    pbaspect([1 1 1])
    ax(2) = subplot(1,2,2);
    imshow(imagesToVisualize.Registered)
    colormap(ax(2), cmap)
    c = colorbar;
    c.Ticks = Ticks/max(Ticks);
    c.TickLabels = TicksLabels;
    ylabel(c,"Trials");
    title("Registered Images Composite Image")
    pbaspect([1 1 1])

    exportgraphics(f,selpath + "/1. Registered Composite.jpg",'Resolution',400)
    close(f)
            
    % Export Registration Coefficients
    f = figure('units','normalized','outerposition',[0 0 1 1], 'Visible','off');
    sgtitle("Pixels Registered")
    s1 = subplot(211);
    plot(regCoeffs(2,:))
    s2 = subplot(212);
    plot(regCoeffs(1,:))
    grid([s1, s2], "minor")
    xlabel([s1, s2], "Frames")
    ylabel([s1, s2], "Pixels")
    title(s1, "Y-Axis Registration")
    title(s2,"X-Axis Registration")

    exportgraphics(f,selpath + "/2. Frames Shifting.jpg",'Resolution',400)
    close(f)

    % All Synapses
    f = figure('units','normalized','outerposition',[0 0 1 1], 'Visible','off');
    M = rgb2gray(imagesToVisualize.Registered);
    bkgArea = M(vertices(1,2):vertices(2,2),vertices(1,1):vertices(3,1));
    bkgMean = mean(bkgArea(:));
    imshow(M,[bkgMean, max(imagesToVisualize.Registered(:))])
    hold on
    gscatter(SynapsesCoG(:,2), SynapsesCoG(:,3), SynapsesCoG(:,6), cell2mat(trilaColors), 'o', 2);
    xlim([0 size(deltaFStacks{1},1)])
    ylim([0 size(deltaFStacks{1},2)])
    grid("minor")
    legend('Selected Synapses', 'Discarded Synapses', 'Location', 'best', 'FontSize',2)
    title("Selected and Discarded Synapses")
    set(gca, 'YDir','reverse')

    exportgraphics(f,selpath + "/3. Selected and Discarded Synapses.jpg",'Resolution',400)
    close(f)

    % Selected Synapses
    nrTrialsTriggered = zeros(length(selectedSynapses.CoG), 1);
    for i = 1 : length(selectedSynapses.CoG)
        nrTrialsTriggered(i) = length(selectedSynapses.trials{i,2});
    end
    f = figure('units','normalized','outerposition',[0 0 1 1], 'Visible','off');
    imshow(M,[bkgMean, max(imagesToVisualize.Registered,[],"all")])
    hold on
    gscatter(selectedSynapses.CoG(:,2), selectedSynapses.CoG(:,3),...
        nrTrialsTriggered*(100/stacksNumber), jet(size(unique(nrTrialsTriggered),1)), "o", 2);
    xlim([0 size(deltaFStacks{1},1)])
    ylim([0 size(deltaFStacks{1},2)])
    grid("minor")
    leg = legend(string(unique(nrTrialsTriggered)*(100/stacksNumber))+"%", 'Location', 'best', 'FontSize',2);
    title("Selected Synapses")
    title(leg, 'Release Probability')
    set(gca, 'YDir','reverse')

    exportgraphics(f,selpath + "/4. Selected Synapses.jpg",'Resolution',400)
    close(f)

    % Synapse Areas
    f = figure('units','normalized','outerposition',[0 0 1 1], 'Visible','off');
    histogram(selectedSynapses.CoG(:,6), numel(unique(selectedSynapses.CoG(:,6))))
    grid("minor")
    xlabel('Area [\mum]','interpreter','Tex')
    ylabel("Counts")
    title("Synapses Areas")

    exportgraphics(f,selpath + "/5. Synapse Areas Histogram.jpg",'Resolution',400)
    close(f)

    % Release Probabilities
    f = figure('units','normalized','outerposition',[0 0 1 1], 'Visible','off');
    hold on
    nrTrialsTriggered = nrTrialsTriggered/stacksNumber*100;
    labels = (1:stacksNumber) / stacksNumber * 100;
    info = categorical(nrTrialsTriggered(nrTrialsTriggered == (100/stacksNumber)),(100/stacksNumber),"Ordinal",true);
    histogram(info)
    info = categorical(nrTrialsTriggered(nrTrialsTriggered > (100/stacksNumber)),labels,"Ordinal",true);
    histogram(info)
    xlabel("(%) Release Probability"); ylabel("Counts")
    h = legend("Responses in a single trial", "Responses in multiple trials");
    pos = h.Position;
    pos(2) = pos(2) - .15;
    annotation('textbox',pos,'String',"Mean Release Probability: " +...
        string(mean(nrTrialsTriggered)) + "%",'FitBoxToText','on', 'Color','r', 'FontSize', 8);
    grid on; grid minor; axis tight

    exportgraphics(f,selpath + "/6. Release Probabilities Histogram.jpg",'Resolution',400)
    close(f)

    % Synapse Traces and Amplitude
    % Only Responses
    f = figure('units','normalized','outerposition',[0 0 1 1], 'Visible','off'); hold on;
    allSynapses = []; 
    for t = 1 : length(unbleachedNormalizedROIs)
        allSynapses(end+1,:) = unbleachedNormalizedROIs{t}(unbleachedNormalizedROIs{t}(:,1) == 3, 2:end);
    end
    avg = mean(allSynapses);
    M = max(avg);
    plot(allSynapses', 'Color', trilaColors{1});
    yline(M, 'g');
    annotation('textbox',[.25 .24 .2 .2], 'String',...
        "Amplitude = " + M, 'FitBoxToText', 'on');
    axis tight; grid minor; grid on; xlabel("time"); ylabel("Amplitude");
    title(["All Synapse Traces", "Only Response Trials Average"])
    exportgraphics(f,selpath + "/7. Synapse Traces and average.jpg",'Resolution',400)
    close(f)
    % Including Failures
    f = figure('units','normalized','outerposition',[0 0 1 1], 'Visible','off'); hold on;
    allSynapses = zeros(length(unbleachedNormalizedROIs), size(unbleachedNormalizedROIs{1},2)-1); 
    for t = 1 : length(unbleachedNormalizedROIs)
        allSynapses(t,:) = unbleachedNormalizedROIs{t}(unbleachedNormalizedROIs{t}(:,1) == 2, 2:end);
    end
    avg = mean(allSynapses);
    M = max(avg);
    plot(allSynapses', 'Color', trilaColors{1});
    yline(M, 'g');
    annotation('textbox',[.25 .24 .2 .2], 'String',...
        "Amplitude = " + M, 'FitBoxToText', 'on');
    axis tight; grid minor; grid on; xlabel("time"); ylabel("Amplitude");
    title(["All Synapse Traces", "All Trials Average"])
    exportgraphics(f,selpath + "/8. Synapse Traces and average (Failures incl.).jpg",'Resolution',400)
    close(f)

    % Synapse Traces
    synapseNrs = selectedSynapses.CoG(:,1);
    selpath = selpath + "Synapse Traces\";
    mkdir(selpath)
    [row, col] = calculateGridFromStacksNumber(stacksNumber, "SynapsesPlots");
    for R = 1 : length(unbleachedNormalizedROIs)
        f = figure('units','normalized','outerposition',[0 0 1 1], 'Visible','off');
        sgtitle("Synapse nr. " + synapseNrs(R))
        for t = 3 : size(unbleachedNormalizedROIs{R}, 1)
            subplot(row, col, t-2)
            if unbleachedNormalizedROIs{R}(t,1) == 1
                plot(unbleachedNormalizedROIs{R}(t,2:end), 'Color', [0.85,0.325,0.098])
            else
                plot(unbleachedNormalizedROIs{R}(t,2:end), 'Color', [0,0.447,0.741])
            end
            axis tight; grid minor; ylim([0.9, 1.4])
            title("Trial " + string(t-2))
        end
        subplot(row, col, t-1)
        plot(unbleachedNormalizedROIs{R}(1,2:end), 'Color', 'k', 'LineWidth',1)
        axis tight; grid minor; ylim([0.9, 1.4]); title("All Trials Average")
        subplot(row, col, t)
        plot(unbleachedNormalizedROIs{R}(2,2:end), 'Color', [0.85,0.325,0.098], 'LineWidth',1)
        axis tight; grid minor; ylim([0.9, 1.4]); title("Responses only Average")
        subplot(row, col, t+1)
        plot(mean(unbleachedNormalizedROIs{R}(unbleachedNormalizedROIs{R}(:,1) == 0, 2:end),1), 'Color', [0,0.447,0.741], 'LineWidth',1)
        axis tight; grid minor; ylim([0.9, 1.4]); title("Failures only Average")

        exportgraphics(f,selpath + "Synapse nr. " + string(synapseNrs(R)) +".jpg",'Resolution',80)
        close(f)
    end
end

function exportPerformanceData(memConsumption, timeConsumption)
    format longG
    filename = "data\TestFiles\Performance.txt";
    fid = fopen(filename,'wt+');
    fprintf(fid, '%s\n', "1. Images Loading:");
    fprintf(fid, '%s', "CPU Time: ");
    fprintf(fid, '%12.6f\n', timeConsumption(1));
    fprintf(fid, '%s', "Memory Consumption: ");
    fprintf(fid, '%12.2f\n', memConsumption(1));
    fprintf(fid, '%s\n', "");
    fprintf(fid, '%s\n', "2. Images Registration:");
    fprintf(fid, '%s', "CPU Time: ");
    fprintf(fid, '%12.6f\n', timeConsumption(2));
    fprintf(fid, '%s', "Memory Consumption: ");
    fprintf(fid, '%12.2f\n', memConsumption(2));
    fprintf(fid, '%s\n', "");
    fprintf(fid, '%s\n', "3. Images Unbleaching:");
    fprintf(fid, '%s', "CPU Time: ");
    fprintf(fid, '%12.6f\n', timeConsumption(3));
    fprintf(fid, '%s', "Memory Consumption: ");
    fprintf(fid, '%12.2f\n', memConsumption(3));
    fprintf(fid, '%s\n', "");
    fprintf(fid, '%s\n', "4. DF Images Calculation:");
    fprintf(fid, '%s', "CPU Time: ");
    fprintf(fid, '%12.6f\n', timeConsumption(4));
    fprintf(fid, '%s', "Memory Consumption: ");
    fprintf(fid, '%12.2f\n', memConsumption(4));
    fprintf(fid, '%s\n', "");
    fprintf(fid, '%s\n', "5. Images Deconvolution:");
    fprintf(fid, '%s', "CPU Time: ");
    fprintf(fid, '%12.6f\n', timeConsumption(5));
    fprintf(fid, '%s', "Memory Consumption: ");
    fprintf(fid, '%12.2f\n', memConsumption(5));
    fprintf(fid, '%s\n', "");
    fprintf(fid, '%s\n', "6. Potential Events Extraction:");
    fprintf(fid, '%s', "CPU Time: ");
    fprintf(fid, '%12.6f\n', timeConsumption(6));
    fprintf(fid, '%s', "Memory Consumption: ");
    fprintf(fid, '%12.2f\n', memConsumption(6));
    fprintf(fid, '%s\n', "");
    fprintf(fid, '%s\n', "7. Synapses Extraction:");
    fprintf(fid, '%s', "CPU Time: ");
    fprintf(fid, '%12.6f\n', timeConsumption(7));
    fprintf(fid, '%s', "Memory Consumption: ");
    fprintf(fid, '%12.2f\n', memConsumption(7));
    fprintf(fid, '%s\n', "");
    fprintf(fid, '%s\n', "8. Results Export:");
    fprintf(fid, '%s', "CPU Time: ");
    fprintf(fid, '%12.6f\n', timeConsumption(8));
    fprintf(fid, '%s', "Memory Consumption: ");
    fprintf(fid, '%12.2f\n', memConsumption(8));
end