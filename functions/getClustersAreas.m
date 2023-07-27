function [clusterAreas, clustersLabels] = getClustersAreas(clusters, pixelFactorNm, edgeLimitNm, clusteringCase)

    limArea = (edgeLimitNm/pixelFactorNm)^2;
    
    switch clusteringCase
        case 'Release Clusters'
            clusters = cat(1, clusters{:});

            clNumbers = unique(clusters(:,1));
            [clustersLabels, clusterAreas] = deal(zeros(numel(clNumbers),1));
            for rc = 1 : numel(clNumbers)
                clus = clusters(clusters(:, 1) == clNumbers(rc),:);
                try
                    [~,av] = convhull(clus(:, 2:3));
                    clusterAreas(rc) = pixelFactorNm * av;
                    if av > limArea
                        clNumbers(rc)
                        clustersLabels(rc) = 1; 
                    end
                catch
                    dist = pdist(clus(:,2:3));
                    if dist == 1
                        clusterAreas(rc) = pixelFactorNm * 0.2;
                    else 
                        clusterAreas(rc) = pixelFactorNm * 0.1;
                    end
                end
            end
        case 'Synapse Clusters'
            clNumbers = unique(clusters(:,1));
            [clustersLabels, clusterAreas] = deal(zeros(numel(clNumbers),1));

            for rc = 1 : numel(clNumbers)
               clus = clusters(clusters(:, 1) == clNumbers(rc),:);
               try
                    [~,av] = convhull(clus(:, 2:3));
                    clusterAreas(rc) = pixelFactorNm * av;
                    if av > limArea
                        clNumbers(rc)
                        clustersLabels(rc) = 1; 
                    end
               catch
                    dist = pdist(clus(:,2:3));
                    if dist == 1
                        clusterAreas(rc) = pixelFactorNm * 0.2;
                    else 
                        clusterAreas(rc) = pixelFactorNm * 0.1;
                    end
               end
            end
    end
end