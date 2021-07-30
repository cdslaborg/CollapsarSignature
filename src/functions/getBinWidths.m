function binWidths = getBinWidths(binEdges)
    % Function that takes as input an array of logarithmically-spaced
    % histogram bin edges and finds their widths.
    
    nbins = length(binEdges)-1;
    binWidths = zeros(1,nbins);
    for i = 1:nbins
        binWidths(i) = binEdges(i+1) - binEdges(i);
    end
end