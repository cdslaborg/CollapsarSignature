function binReduction = binReductionStairPlot(n,binEdges,lineColor,lineWidth)
    % Function that combines adjacent bins with < 5 events
    
    if nargin < 4; lineWidth = 0.5; end
    i = 1;
    while i < length(n)+1
        if n(i) < 5
            if i == length(n)
                n(i-1) = n(i-1) + n(i);
                n(i) = [];
                binEdges(i) = [];
            else
                n(i) = n(i) + n(i+1);
                n(i+1) = [];
                binEdges(i+1) = [];
                continue
            end
        else
            i = i+1;
        end
    end
    nbins = length(binEdges)-1;
    binCenters = zeros(1,nbins);
    for i = 1:nbins
        binCenters(i) = 0.5*(binEdges(i) + binEdges(i+1));
    end
    dNdT = n ./ binCenters;
    
    % Plot new dN/dT
    stairs(binCenters, dNdT, 'color', lineColor, 'lineWidth', lineWidth);
    
    binReduction = nbins;
end