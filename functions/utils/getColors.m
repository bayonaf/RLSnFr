function [colors, marksSize] = getColors(clustersAverage, palleteSize)
    pallete = getPallete(palleteSize);
    N = size(clustersAverage,1);
    colors = zeros(N, 3);
    marksSize = zeros(N, 1);
    for i = 1 : N
        [x,y] = deal(fix(clustersAverage(i,2)), fix(clustersAverage(i,1)));
        colors(i,:) = pallete(x,y,:);
        X = [fix(palleteSize/2) fix(palleteSize/2); x y];
        marksSize(i,:) = (abs(pdist(X,'euclidean'))+10)/3;
    end
end

