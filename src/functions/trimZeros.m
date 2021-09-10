function [counts, binWidths, binEdges, binCenters] = trimZeros(counts, binEdges)
    % Function that takes as input histogram-like data and trims off all
    % values on the tails of the distribution that are either zeros or are 
    % surrounded by zeros. Assumes zeros are only on the tails.
    
    lenCounts = length(counts);
    if lenCounts ~= length(binEdges)-1; error("'binEdges' is not 1 + length of 'counts'"); end
    
    half = round(lenCounts/2);
    zeroIndicies = find(~counts);
    if ~isempty(zeroIndicies)
        zeroIndiciesLeft = find(zeroIndicies < half);
        zeroIndiciesRight = find(zeroIndicies > half);
        if ~isempty(zeroIndiciesRight)
            zeroIndiciesRight = zeroIndicies(zeroIndiciesRight);
            counts(min(zeroIndiciesRight):end) = [];
            binEdges(min(zeroIndiciesRight)+1:end) = [];
        end
        if ~isempty(zeroIndiciesLeft)
            zeroIndiciesLeft = zeroIndicies(zeroIndiciesLeft);
            counts(1:max(zeroIndiciesLeft)) = [];
            binEdges(1:max(zeroIndiciesLeft)) = [];
        end
    end
    binWidths = getBinWidths(binEdges);
    binCenters = getBinCenters(binEdges);
end