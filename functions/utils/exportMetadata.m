function exportMetadata(selpath, experimentName, threshFactor, threshold, bkgrVertices, bkgrAvgs, xbegin, paramsMap)
% EXPORTMETADATA(): This function exports a txt file containing app
% configuration parameters in execution time defined by the user
%   
%   EXPORTMETADATA(selpath, experimentName, threshFactor, threshold, bkgrVertices, bkgrAvgs, xbegin, paramsMap) 
%   
%   SELPATH:        String containing the exportation URL 
%   EXPERIMENTNAME: String containing the experiment name
%   THRESHFACTOR:   Threshold factor defined by user to get potential
%   release events
%   THRESHOLD:      Numerical Threshold value defined by threshFactor
%   BKGRVERTICES:   4x4 matrix containing the square vertices of the
%   selected background area
%   BKGRAVGS:       Average background values substracted per trial
%   XBEGIN:         Number of initial frames with unusual brightness
%   PARAMSMAP:      Parameters map containing the app config values in the
%   file key_value.m

    filename = selpath + "Metadata.txt";
    
    fid = fopen(filename,'wt+');

    fprintf(fid, '%s\n', "Experiment Name:");
    fprintf(fid, '%s\n', experimentName);

    fprintf(fid, '%s\n', "Threshold Factor:");
    fprintf(fid, '%s\n', threshFactor);

    fprintf(fid, '%s\n', "Threshold Value:");
    fprintf(fid, '%s\n', threshold);

    if ~isempty(bkgrVertices)
        fprintf(fid, '%s\n', "Background Area Vertices [X, Y]:");
        fprintf(fid, '%d %d\n', bkgrVertices.');
    
        fprintf(fid, '%s\n', "Background Substracted Stats:");
        fprintf(fid, '%s\n', "Averaged Background Substracted:");
        fprintf(fid, '%s\n', mean(bkgrAvgs));
        fprintf(fid, '%s\n', "Minimal Background Substracted:");
        fprintf(fid, '%s\n', min(bkgrAvgs));
        fprintf(fid, '%s\n', "Maximal Background Substracted:");
        fprintf(fid, '%s\n', max(bkgrAvgs));
    else
        fprintf(fid, '%s\n', "Background Substracted Defined by User:");
        fprintf(fid, '%s\n', bkgrAvgs);
    end

    fprintf(fid, '%s\n', '');
    fprintf(fid, '%s\n', "EXPERIMENT PARAMETERS");
    fprintf(fid, '%s ', "Signal Trigger Frame:");
    fprintf(fid, '%s\n', string(paramsMap('triggerTime')));
    fprintf(fid, '%s ', "Baseline Initial Frame:");
    fprintf(fid, '%s\n', string(paramsMap('initFrames')));
    fprintf(fid, '%s ', "Baseline Ending Frame:");
    fprintf(fid, '%s\n', string(paramsMap('baselineEndFrame')));
    fprintf(fid, '%s ', "Unusuall Brightness Initial Frame:");
    fprintf(fid, '%s\n', string(xbegin));

    fprintf(fid, '%s\n', "Results Exported to:");
    fprintf(fid, '%s\n', selpath);

    fclose(fid);
end

