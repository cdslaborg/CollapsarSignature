function hist2stairs(counts, binEdges, lineColor, lineWidth)
    % Correctly renders histogram data as a stairs function.
    
    counts(end+1) = counts(end);
    stairs(binEdges, counts, 'color', lineColor, 'lineWidth', lineWidth);
end