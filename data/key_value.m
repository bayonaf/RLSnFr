% EXPERIMENTTYPE TYPE (Spontaneous = 1 or Triggered = 0)
experimentType=0
% Variables for RL Deconvolution
% function startupFcn(app) Deconv_Sniffer
RLIterations=1000
triggerTime=21
% Variables for Release/synapse Clusters (Distances in nanometers)
% Function getClusterIndexes()
epsilonRCnm=246
epsilonSCnm=246
% Proximity radius for Synapse Clusters is 4 times the gaussian sigma above
radiusSCnm=1600
pixelFactorNm=246
minPts=1
% Thresholding factor
% function startupFcn(app) Results_Analyzis
threshFactor=3.5
% Initial frames with unexpected brigthness
initFrames=3
baselineEndFrame=17
% Used in Gaussian averaging
gaussSigmaSA=400
% Used for signals pseudofitting
pseudoTriggerTime=10
% Max distance from branches area for sigmoid balance
maxBranchesDistNm=1000
% Registration Mode 0='monomodal', 1 = 'multimodal'
regMode=0
% FOR SPONTANEOUS RELASE DECONVOLUTION
numStacks=20
windowOffset=10
windowLength=15
overlappedFrames=10
stdFactor=4