function [lower, upper] = plotLimits(data, type, optionalData)
    % Determines the axis limits for a plot of 'data' with a 10% margin.
    % Can optionally exclude values of 'data' where 'optionalData' is zero.
    
    if ~any(strcmpi(type,["linear", "log"])); error("'type' is not an allowable value"); end
    margin = 0.1;
    
    if nargin == 2
        minimum = min(data);
        maximum = max(data);
        if strcmpi(type, "linear")
            lower = minimum - (maximum - minimum) * margin;
            upper = maximum + (maximum - minimum) * margin;
        else
            logmin = log(min(data(data > 0)));
            logmax = log(maximum);
            lower = exp(logmin - (logmax - logmin) * margin);
            upper = exp(logmax + (logmax - logmin) * margin);
        end
        
    elseif nargin == 3
        if length(data) ~= length(optionalData); error("data sets are not the same size"); end
        minimum = min(data(optionalData > 0));
        maximum = max(data(optionalData > 0));
        if strcmpi(type, "linear")
            lower = minimum - (maximum - minimum) * margin;
            upper = maximum + (maximum - minimum) * margin;
        else
            logmin = log(min(data(data > 0 & optionalData > 0)));
            logmax = log(maximum);
            lower = exp(logmin - (logmax - logmin) * margin);
            upper = exp(logmax + (logmax - logmin) * margin);
        end
    end
 end