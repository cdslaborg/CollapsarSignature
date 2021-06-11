function binCenters = getBinCenters(binEdges)
    % Function that takes as input an array of logarithmically-spaced
    % histogram bin edges and finds their centers in log space.
    
    nbins = length(binEdges)-1;
    binCenters = zeros(1,nbins);
    for i = 1:nbins
        binCenters(i) = exp(0.5*(log(binEdges(i)) + log(binEdges(i+1))));
    end
end