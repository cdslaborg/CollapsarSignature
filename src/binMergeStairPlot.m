function nbins = binMergeStairPlot(counts,binEdges,lineColor,lineWidth)
    % Function that combines adjacent bins with < 5 events
    
    if nargin < 4; error("Not enough input arguments"); end
    
    % Merge bins from left to right
    half = round(length(counts)/2);
    i = 1;
    while i < half
        if counts(i) < 5
            counts(i) = counts(i) + counts(i+1);
            counts(i+1) = [];
            binEdges(i+1) = [];
            half = half - 1;
            continue
        else
            i = i+1;
        end
    end
    
    % Merge bins from right to left
    i = 0;
    while i < length(counts)-half
        if counts(end-i) < 5
            counts(end-1-i) = counts(end-1-i) + counts(end-i);
            counts(end-i) = [];
            binEdges(end-i-1) = [];
        else
            i = i+1;
        end
    end
    
    nbins = length(binEdges)-1;
    binCenters = getBinCenters(binEdges);
    dNdT = counts ./ binCenters;
    
    % Plot new dN/dT
    hist2stairs(dNdT, binEdges, lineColor, lineWidth)
end