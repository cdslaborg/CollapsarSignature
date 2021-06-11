function binEdges = getBinEdges(data,nbins)
    % Function that generates the bin edges for 'nbins' bins that are
    % equally spaced in log space with limits 1% below and 1% above 'data'.
    
    binEdges = exp(linspace(log(min(data)*99/100),log(max(data)*101/100),nbins+1));
end