clear all;
close all;
filePath = mfilename('fullpath');
[currentDir,fileName,fileExt] = fileparts(filePath); cd(currentDir);
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath([genpath("../../../libmatlab"), genpath("./functions")], "-begin");

% Figure Parameters
fontSize = 14;
lineWidth = 1.5;
figureColor = "white";
defaultColor = [0, 0.4470, 0.7410];
secondaryColor = '#FFA500'; %orange
sgrbColor = "blue";
lgrbColor = "red";
allColor = "magenta";
confidenceColor = [0.5, 0.5, 0.5]; %gray

% Simulation Options
grbType = "all"; % choose between "lgrb", "sgrb", and "all"
freshRunEnabled = false; % this must be set to 'true' for first ever simulation. Thereafter, it can be set to 'false' to save time.
saveNewImages = false; % export figure on or off. Must have export_fig installed.
skip.lgrb = 28; % reduce synthetic lgrb dataset of 200,000 to ~1400, equivalent to BATSE; 'freshRunEnabled' must be 'true'
skip.sgrb = 10; % reduce synthetic sgrb dataset of 18,848 to ~600, equivalent to BATSE; 'freshRunEnabled' must be 'true'
numBins = 59; % number of bins to bin synthetic data in; 'freshRunEnabled' must be 'true'
numRuns = 10000; % number of detector simulation runs; 'freshRunEnabled' must be 'true'
matFileName = "synSamT90.mat";

if ~any(strcmpi(grbType,["lgrb","sgrb","all"])); error("'grbType' is not an allowable value"); end
if freshRunEnabled
    % Begin Data Processing
    lgrb = importdata("..\in\syntheticSampleB10.csv");
    sgrb = importdata("..\in\SyntheticSample_All.txt");
    
    icol = struct();
    icol.lgrb.logT90 = 8; % column index of logDur
    icol.lgrb.detProb = 10; % column index of detection probability
    icol.sgrb.T90 = 9; % column index of T90
    icol.sgrb.detProb = 23; % column index of detection probability
    
    synSam = struct(); 
    synSam.lgrb.T90 = exp(lgrb.data(:,icol.lgrb.logT90));
    synSam.lgrb.detProb = lgrb.data(:,icol.lgrb.detProb);
    synSam.sgrb.T90 = sgrb.data(:,icol.sgrb.T90);
    synSam.sgrb.detProb = sgrb.data(:,icol.sgrb.detProb);
    synSam.all.binEdges = getBinEdges([synSam.sgrb.T90;synSam.lgrb.T90],numBins);
    synSam.all.binCenters = getBinCenters(synSam.all.binEdges);

    synSam.lgrb.whole.ndata = length(synSam.lgrb.T90);
    synSam.lgrb.whole.countsVec = histcounts(synSam.lgrb.T90,synSam.all.binEdges);
    [synSam.lgrb.whole.countsVec, synSam.lgrb.whole.binCenters, synSam.lgrb.whole.binEdges] = trimZeros(synSam.lgrb.whole.countsVec, synSam.all.binEdges);
    synSam.sgrb.whole.ndata = length(synSam.sgrb.T90);
    synSam.sgrb.whole.countsVec = histcounts(synSam.sgrb.T90,synSam.all.binEdges);
    [synSam.sgrb.whole.countsVec, synSam.sgrb.whole.binCenters, synSam.sgrb.whole.binEdges] = trimZeros(synSam.sgrb.whole.countsVec, synSam.all.binEdges);
    synSam.all.whole.ndata = length([synSam.lgrb.T90;synSam.sgrb.T90]);
    synSam.all.whole.countsVec = histcounts([synSam.lgrb.T90;synSam.sgrb.T90],synSam.all.binEdges);    
    
    % Normalize the ratio of LGRBs/SGRBs to the observed ratio for BATSE in the combined, pre-detected dataset. Does not affect any other computations.
    synSam.all.whole.normalized.ndata = length([synSam.lgrb.T90(1:round(length(synSam.sgrb.T90)*1366/600));synSam.sgrb.T90]);
    synSam.all.whole.normalized.countsVec = histcounts([synSam.lgrb.T90(1:round(length(synSam.sgrb.T90)*1366/600));synSam.sgrb.T90],synSam.all.binEdges);
    [synSam.all.whole.normalized.countsVec, synSam.all.whole.normalized.binCenters, synSam.all.whole.normalized.binEdges] ...
        = trimZeros(synSam.all.whole.normalized.countsVec, synSam.all.binEdges);
    
    % If detection probability is above a uniformly generated random number, observe event. Repeat 'numRuns' times.
    for i = 1:numRuns
        Mask.lgrb = lgrb.data(:,icol.lgrb.detProb) > unifrnd(0,1,length(lgrb.data),1);
        Mask.sgrb = sgrb.data(:,icol.sgrb.detProb) > unifrnd(0,1,length(sgrb.data),1);
        T90.lgrb = exp(lgrb.data(Mask.lgrb,icol.lgrb.logT90));  % choose those events that would be detected
        T90.sgrb = sgrb.data(Mask.sgrb,icol.sgrb.T90);          % choose those events that would be detected
        T90.lgrb = T90.lgrb(1:skip.lgrb:end);   % reduce the detected dataset
        T90.sgrb = T90.sgrb(1:skip.sgrb:end);   % reduce the detected dataset
        synSam.lgrb.reducedVec(i).ndata = length(T90.lgrb);
        synSam.lgrb.reducedVec(i).countsVec = histcounts(T90.lgrb,synSam.all.binEdges);
        synSam.sgrb.reducedVec(i).ndata = length(T90.sgrb);
        synSam.sgrb.reducedVec(i).countsVec = histcounts(T90.sgrb,synSam.all.binEdges);
        synSam.all.reducedVec(i).ndata = length([T90.lgrb;T90.sgrb]);
        synSam.all.reducedVec(i).countsVec = histcounts([T90.lgrb;T90.sgrb],synSam.all.binEdges);
    end
    clear('lgrb','sgrb','icol','Mask','T90');
    
    for i = 1:numBins
        % Switch the positions of countsVec and reducedVec in the struct
        for j = 1:numRuns
            synSam.lgrb.countsVec(i).reducedVec(j) = synSam.lgrb.reducedVec(j).countsVec(i);
            synSam.sgrb.countsVec(i).reducedVec(j) = synSam.sgrb.reducedVec(j).countsVec(i);
            synSam.all.countsVec(i).reducedVec(j) = synSam.all.reducedVec(j).countsVec(i);
        end
        % Calculate the 95th percentile, the mean, and the 5th percentile
        synSam.lgrb.reduced.percentile95.countsVec(i) = round(prctile(synSam.lgrb.countsVec(i).reducedVec,95));
        synSam.lgrb.reduced.percentile50.countsVec(i) = round(prctile(synSam.lgrb.countsVec(i).reducedVec,50));
        synSam.lgrb.reduced.percentile05.countsVec(i) = round(prctile(synSam.lgrb.countsVec(i).reducedVec,5));
        synSam.sgrb.reduced.percentile95.countsVec(i) = round(prctile(synSam.sgrb.countsVec(i).reducedVec,95));
        synSam.sgrb.reduced.percentile50.countsVec(i) = round(prctile(synSam.sgrb.countsVec(i).reducedVec,50));
        synSam.sgrb.reduced.percentile05.countsVec(i) = round(prctile(synSam.sgrb.countsVec(i).reducedVec,5));
        synSam.all.reduced.percentile95.countsVec(i) = round(prctile(synSam.all.countsVec(i).reducedVec,95));
        synSam.all.reduced.percentile50.countsVec(i) = round(prctile(synSam.all.countsVec(i).reducedVec,50));
        synSam.all.reduced.percentile05.countsVec(i) = round(prctile(synSam.all.countsVec(i).reducedVec,5));
    end
    clear('i','j');
    synSam.lgrb = rmfield(synSam.lgrb, {'countsVec', 'reducedVec'});
    synSam.sgrb = rmfield(synSam.sgrb, {'countsVec', 'reducedVec'});
    synSam.all = rmfield(synSam.all, {'countsVec', 'reducedVec'});
    
    % Calculate the mean value of ndata across all trials
    synSam.lgrb.reduced.percentile50.ndata = sum(synSam.lgrb.reduced.percentile50.countsVec);
    synSam.sgrb.reduced.percentile50.ndata = sum(synSam.sgrb.reduced.percentile50.countsVec);
    synSam.all.reduced.percentile50.ndata = sum(synSam.all.reduced.percentile50.countsVec);
    
    % Trim zeros off of binned data for later curve fitting
    [synSam.lgrb.reduced.percentile95.countsVec, synSam.lgrb.reduced.percentile95.binCenters, synSam.lgrb.reduced.percentile95.binEdges] ...
        = trimZeros(synSam.lgrb.reduced.percentile95.countsVec, synSam.all.binEdges);
    [synSam.lgrb.reduced.percentile50.countsVec, synSam.lgrb.reduced.percentile50.binCenters, synSam.lgrb.reduced.percentile50.binEdges] ...
        = trimZeros(synSam.lgrb.reduced.percentile50.countsVec, synSam.all.binEdges);
    [synSam.lgrb.reduced.percentile05.countsVec, synSam.lgrb.reduced.percentile05.binCenters, synSam.lgrb.reduced.percentile05.binEdges] ...
        = trimZeros(synSam.lgrb.reduced.percentile05.countsVec, synSam.all.binEdges);
    [synSam.sgrb.reduced.percentile95.countsVec, synSam.sgrb.reduced.percentile95.binCenters, synSam.sgrb.reduced.percentile95.binEdges] ...
        = trimZeros(synSam.sgrb.reduced.percentile95.countsVec, synSam.all.binEdges);
    [synSam.sgrb.reduced.percentile50.countsVec, synSam.sgrb.reduced.percentile50.binCenters, synSam.sgrb.reduced.percentile50.binEdges] ...
        = trimZeros(synSam.sgrb.reduced.percentile50.countsVec, synSam.all.binEdges);
    [synSam.sgrb.reduced.percentile05.countsVec, synSam.sgrb.reduced.percentile05.binCenters, synSam.sgrb.reduced.percentile05.binEdges] ...
        = trimZeros(synSam.sgrb.reduced.percentile05.countsVec, synSam.all.binEdges);
    [synSam.all.reduced.percentile95.countsVec, synSam.all.reduced.percentile95.binCenters, synSam.all.reduced.percentile95.binEdges] ...
        = trimZeros(synSam.all.reduced.percentile95.countsVec, synSam.all.binEdges);
    [synSam.all.reduced.percentile50.countsVec, synSam.all.reduced.percentile50.binCenters, synSam.all.reduced.percentile50.binEdges] ...
        = trimZeros(synSam.all.reduced.percentile50.countsVec, synSam.all.binEdges);
    [synSam.all.reduced.percentile05.countsVec, synSam.all.reduced.percentile05.binCenters, synSam.all.reduced.percentile05.binEdges] ...
        = trimZeros(synSam.all.reduced.percentile05.countsVec, synSam.all.binEdges);
    
    % Convert the dependent variable dN/dlog(T) to dN/dT by dividing it by T90
    synSam.lgrb.whole.dNdT = synSam.lgrb.whole.countsVec ./ synSam.lgrb.whole.binCenters;
    synSam.sgrb.whole.dNdT = synSam.sgrb.whole.countsVec ./ synSam.sgrb.whole.binCenters;
    synSam.all.whole.dNdT = synSam.all.whole.countsVec ./ synSam.all.binCenters;
    synSam.all.whole.normalized.dNdT = synSam.all.whole.normalized.countsVec ./ synSam.all.whole.normalized.binCenters;
    synSam.lgrb.reduced.percentile50.dNdT = synSam.lgrb.reduced.percentile50.countsVec ./ synSam.lgrb.reduced.percentile50.binCenters;
    synSam.sgrb.reduced.percentile50.dNdT = synSam.sgrb.reduced.percentile50.countsVec ./ synSam.sgrb.reduced.percentile50.binCenters;
    synSam.all.reduced.percentile50.dNdT = synSam.all.reduced.percentile50.countsVec ./ synSam.all.reduced.percentile50.binCenters;
    
    % Create a smooth fit for the binned data. Must have Curve Fitting Toolbox installed.
    fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,1,1]);
    ft = fittype('a*exp(-((log(x)-b)/c)^2)', 'options', fo); % fit a Gaussian
    synSam.lgrb.reduced.percentile50.f = fit(synSam.lgrb.reduced.percentile50.binCenters.', synSam.lgrb.reduced.percentile50.countsVec.', ft);
    synSam.lgrb.reduced.percentile50.fdNdT = synSam.lgrb.reduced.percentile50.f(synSam.lgrb.reduced.percentile50.binCenters).' ...
        ./ synSam.lgrb.reduced.percentile50.binCenters;
    synSam.lgrb.reduced.percentile95.f = fit(synSam.lgrb.reduced.percentile95.binCenters.', synSam.lgrb.reduced.percentile95.countsVec.', ft);
    synSam.lgrb.reduced.percentile95.fdNdT = synSam.lgrb.reduced.percentile95.f(synSam.lgrb.reduced.percentile95.binCenters).' ...
        ./ synSam.lgrb.reduced.percentile95.binCenters;
    synSam.lgrb.reduced.percentile05.f = fit(synSam.lgrb.reduced.percentile05.binCenters.', synSam.lgrb.reduced.percentile05.countsVec.', ft);
    synSam.lgrb.reduced.percentile05.fdNdT = synSam.lgrb.reduced.percentile05.f(synSam.lgrb.reduced.percentile05.binCenters).' ...
        ./ synSam.lgrb.reduced.percentile05.binCenters;
    synSam.sgrb.reduced.percentile50.f = fit(synSam.sgrb.reduced.percentile50.binCenters.', synSam.sgrb.reduced.percentile50.countsVec.', ft);
    synSam.sgrb.reduced.percentile50.fdNdT = synSam.sgrb.reduced.percentile50.f(synSam.sgrb.reduced.percentile50.binCenters).' ...
        ./ synSam.sgrb.reduced.percentile50.binCenters;
    synSam.sgrb.reduced.percentile95.f = fit(synSam.sgrb.reduced.percentile95.binCenters.', synSam.sgrb.reduced.percentile95.countsVec.', ft);
    synSam.sgrb.reduced.percentile95.fdNdT = synSam.sgrb.reduced.percentile95.f(synSam.sgrb.reduced.percentile95.binCenters).' ...
        ./ synSam.sgrb.reduced.percentile95.binCenters;
    synSam.sgrb.reduced.percentile05.f = fit(synSam.sgrb.reduced.percentile05.binCenters.', synSam.sgrb.reduced.percentile05.countsVec.', ft);
    synSam.sgrb.reduced.percentile05.fdNdT = synSam.sgrb.reduced.percentile05.f(synSam.sgrb.reduced.percentile05.binCenters).' ...
        ./ synSam.sgrb.reduced.percentile05.binCenters;
    synSam.all.reduced.percentile50.fdNdT = synSam.lgrb.reduced.percentile50.f(synSam.all.reduced.percentile50.binCenters).' ...
        ./ synSam.all.reduced.percentile50.binCenters + synSam.sgrb.reduced.percentile50.f(synSam.all.reduced.percentile50.binCenters).' ...
        ./ synSam.all.reduced.percentile50.binCenters;
    
    % Save mat file
    synSam.output.path = "../out"; if ~isfolder(synSam.output.path); mkdir(synSam.output.path); end
    save(synSam.output.path + "/" + matFileName,"synSam");
    
else
    % Load mat file
    synSam.output.path = "../out";
    load(synSam.output.path + "/" + matFileName); % loads synSam object
end

% Figure 1: Histogram for whole dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        histogram('BinEdges', synSam.lgrb.whole.binEdges, 'BinCounts', synSam.lgrb.whole.countsVec);
        legend(["Whole LGRB synthetic dataset," + newline + "n = " + synSam.lgrb.whole.ndata + ", in " ...
               + length(synSam.lgrb.whole.countsVec) + " bins"], "interpreter", "tex", "location", "northwest", "fontsize", fontSize-3);
        ylim([1, 3e5]);
    elseif strcmpi(grbType,"sgrb")
        histogram('BinEdges',synSam.sgrb.whole.binEdges, 'BinCounts', synSam.sgrb.whole.countsVec);
        legend(["Whole SGRB synthetic dataset," + newline + "n = " + synSam.sgrb.whole.ndata + ", in " ...
               + length(synSam.sgrb.whole.countsVec) + " bins"], "interpreter", "tex", "location", "northwest", "fontsize", fontSize-3);
        ylim([1, 1e4]);
    elseif strcmpi(grbType,"all")
        histogram('BinEdges', synSam.all.whole.normalized.binEdges, 'BinCounts', synSam.all.whole.normalized.countsVec, 'FaceColor', defaultColor);
        %xline(2, 'color', 'green', 'lineWidth', lineWidth);
        legend(["Whole, BATSE-GRB-ratio-normalized" + newline + "synthetic dataset, n = " + synSam.all.whole.normalized.ndata + ", in " ...
               + length(synSam.all.whole.normalized.countsVec) + " bins"], "interpreter", "tex", "location", "northwest", "fontsize", fontSize-3);
        [xlower, xupper] = plotLimits(synSam.all.whole.normalized.binEdges, 'log');
        xlim([xlower, xupper]);
        ylim([1, 1e5]);
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dlog(T_{90})", "interpreter", "tex", "fontsize", fontSize);
    if saveNewImages
        if strcmpi(grbType,"lgrb")
            export_fig(synSam.output.path + "/" + "lgrbHistogramWhole.png", "-m4 -transparent");
        elseif strcmpi(grbType,"sgrb")
            export_fig(synSam.output.path + "/" + "sgrbHistogramWhole.png", "-m4 -transparent");
        elseif strcmpi(grbType,"all")
            export_fig(synSam.output.path + "/" + "allHistogramWhole.png", "-m4 -transparent");
        end
    end
hold off;

% Figure 2: Plot dN/dT of whole dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        hist2stairs(synSam.lgrb.whole.dNdT, synSam.lgrb.whole.binEdges, defaultColor, lineWidth);
        legend(["Whole LGRB synthetic dataset," + newline + "n = " + synSam.lgrb.whole.ndata + ", in " ...
               + length(synSam.lgrb.whole.dNdT) + " bins"], "interpreter", "tex", "location", "southwest", "fontsize", fontSize-3);
    elseif strcmpi(grbType,"sgrb")
        hist2stairs(synSam.sgrb.whole.dNdT, synSam.sgrb.whole.binEdges, defaultColor, lineWidth);
        legend(["Whole SGRB synthetic dataset," + newline + "n = " + synSam.sgrb.whole.ndata + ", in " ...
               + length(synSam.sgrb.whole.dNdT) + " bins"], "interpreter", "tex", "location", "southwest", "fontsize", fontSize-3);
        ylim([1e-2, 1e5]);
    elseif strcmpi(grbType,"all")
        hist2stairs(synSam.all.whole.normalized.dNdT, synSam.all.whole.normalized.binEdges, defaultColor, lineWidth);
        % new
        %nbins = binMergeStairPlot(synSam.all.whole.normalized.countsVec, synSam.all.whole.normalized.binEdges, defaultColor, lineWidth);
        % end new
        legend(["Whole, BATSE-GRB-ratio-normalized" + newline + "synthetic dataset, n = " + synSam.all.whole.normalized.ndata + ", in " ...
               + length(synSam.all.whole.normalized.countsVec) + " bins"], "interpreter", "tex", "location", "southwest", "fontsize", fontSize-3);
        xlim([xlower, xupper]);
        [ylower, yupper] = plotLimits(synSam.all.whole.normalized.dNdT, 'log');
        ylim([ylower, yupper]);
        yticks(10.^(linspace(-4,4,5)));
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
    if saveNewImages
        if strcmpi(grbType,"lgrb")
            export_fig(synSam.output.path + "/" + "lgrbDNDTWhole.png", "-m4 -transparent");
        elseif strcmpi(grbType,"sgrb")
            export_fig(synSam.output.path + "/" + "sgrbDNDTWhole.png", "-m4 -transparent");
        elseif strcmpi(grbType,"all")
            export_fig(synSam.output.path + "/" + "allDNDTWhole.png", "-m4 -transparent");
        end
    end
hold off;

% Figure 3: Plot dN/dT of reduced dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        plot(synSam.lgrb.reduced.percentile50.binCenters, synSam.lgrb.reduced.percentile50.fdNdT, 'color', lgrbColor, 'lineWidth', lineWidth);
        hist2stairs(synSam.lgrb.reduced.percentile50.dNdT, synSam.lgrb.reduced.percentile50.binEdges, lgrbColor, lineWidth);
        plot(synSam.lgrb.reduced.percentile95.binCenters, synSam.lgrb.reduced.percentile95.fdNdT, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        plot(synSam.lgrb.reduced.percentile05.binCenters, synSam.lgrb.reduced.percentile05.fdNdT, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        legend([synSam.lgrb.reduced.percentile50.ndata + " synth. observed LGRBs" + newline + "spanning " + length(synSam.lgrb.reduced.percentile50.fdNdT) ...
               + " bins", "", "90% confidence Interval", ""], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    elseif strcmpi(grbType,"sgrb")
        plot(synSam.sgrb.reduced.percentile50.binCenters, synSam.sgrb.reduced.percentile50.fdNdT, 'color', sgrbColor, 'lineWidth', lineWidth);
        hist2stairs(synSam.sgrb.reduced.percentile50.dNdT, synSam.sgrb.reduced.percentile50.binEdges, sgrbColor, lineWidth);
        plot(synSam.sgrb.reduced.percentile95.binCenters, synSam.sgrb.reduced.percentile95.fdNdT, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        plot(synSam.sgrb.reduced.percentile05.binCenters, synSam.sgrb.reduced.percentile05.fdNdT, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        legend([synSam.sgrb.reduced.percentile50.ndata + " synth. observed SGRBs" + newline + "spanning " + length(synSam.sgrb.reduced.percentile50.fdNdT) ...
               + " bins", "", "90% confidence Interval", ""], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    elseif strcmpi(grbType,"all")
        plot(synSam.all.reduced.percentile50.binCenters, synSam.all.reduced.percentile50.fdNdT*10, 'color', 'magenta', 'lineWidth', lineWidth);
        hist2stairs(synSam.lgrb.reduced.percentile50.dNdT, synSam.lgrb.reduced.percentile50.binEdges, lgrbColor, lineWidth-0.5);
        plot(synSam.lgrb.reduced.percentile50.binCenters, synSam.lgrb.reduced.percentile50.fdNdT, 'color', lgrbColor, 'lineWidth', lineWidth);
        hist2stairs(synSam.sgrb.reduced.percentile50.dNdT, synSam.sgrb.reduced.percentile50.binEdges, sgrbColor, lineWidth-0.5);
        plot(synSam.sgrb.reduced.percentile50.binCenters, synSam.sgrb.reduced.percentile50.fdNdT, 'color', sgrbColor, 'lineWidth', lineWidth);
        plot(synSam.lgrb.reduced.percentile95.binCenters, synSam.lgrb.reduced.percentile95.fdNdT, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        plot(synSam.lgrb.reduced.percentile05.binCenters, synSam.lgrb.reduced.percentile05.fdNdT, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        plot(synSam.sgrb.reduced.percentile95.binCenters, synSam.sgrb.reduced.percentile95.fdNdT, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        plot(synSam.sgrb.reduced.percentile05.binCenters, synSam.sgrb.reduced.percentile05.fdNdT, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        xlim([1e-2, 2e3]);
        %ylim([1e-4, 1e3]);
        legend([synSam.all.reduced.percentile50.ndata + " synth. observed GRBs" + newline + "spanning " + length(synSam.all.reduced.percentile50.fdNdT) + ...
               " bins, x10", "", synSam.lgrb.reduced.percentile50.ndata + " synth. observed LGRBs", "", synSam.sgrb.reduced.percentile50.ndata + ...
               " synth. observed SGRBs", "90% confidence Interval"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        %legend('boxoff');
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
    if saveNewImages
        if strcmpi(grbType,"lgrb")
            export_fig(synSam.output.path + "/" + "lgrbDNDTObserved.png", "-m4 -transparent");
        elseif strcmpi(grbType,"sgrb")
            export_fig(synSam.output.path + "/" + "sgrbDNDTObserved.png", "-m4 -transparent");
        elseif strcmpi(grbType,"all")
            export_fig(synSam.output.path + "/" + "allDNDTObserved.png", "-m4 -transparent");
        end
    end
hold off;

% Figure 4: Combine adjacent bins with fewer than 5 events and plot dN/dT of reduced dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        hist2stairs(synSam.lgrb.reduced.percentile50.dNdT, synSam.lgrb.reduced.percentile50.binEdges, secondaryColor, lineWidth);
        nbins = binMergeStairPlot(synSam.lgrb.reduced.percentile50.countsVec, synSam.lgrb.reduced.percentile50.binEdges, defaultColor, lineWidth);
        legend([synSam.lgrb.reduced.percentile50.ndata + " synth. observed LGRBs" + newline + "spanning " + length(synSam.lgrb.reduced.percentile50.dNdT) ...
               + " bins", "Adjacent bins with <5 events merged," + newline + "resulting in " + nbins + " bins"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    elseif strcmpi(grbType,"sgrb")
        hist2stairs(synSam.sgrb.reduced.percentile50.dNdT, synSam.sgrb.reduced.percentile50.binEdges, secondaryColor, lineWidth);
        nbins = binMergeStairPlot(synSam.sgrb.reduced.percentile50.countsVec, synSam.sgrb.reduced.percentile50.binEdges, defaultColor, lineWidth);
        legend([synSam.sgrb.reduced.percentile50.ndata + " synth. observed SGRBs" + newline + "spanning " + length(synSam.sgrb.reduced.percentile50.dNdT) ...
               + " bins", "Adjacent bins with <5 events merged," + newline + "resulting in " + nbins + " bins"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    elseif strcmpi(grbType,"all")
        hist2stairs(synSam.all.reduced.percentile50.dNdT, synSam.all.reduced.percentile50.binEdges, defaultColor, lineWidth);
        nbins = binMergeStairPlot(synSam.all.reduced.percentile50.countsVec, synSam.all.reduced.percentile50.binEdges, secondaryColor, lineWidth);
        legend([synSam.all.reduced.percentile50.ndata + " synth. observed GRBs" + newline + "spanning " + length(synSam.all.reduced.percentile50.dNdT) + ...
               " bins", "Adjacent bins with <5 events merged," + newline + "resulting in " + nbins + " bins"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        xlim([1e-2, 2e3]);
        ylim([5e-4, 1e3]);
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
    if saveNewImages
        if strcmpi(grbType,"lgrb")
            export_fig(synSam.output.path + "/" + "lgrbDNDTObservedMerged.png", "-m4 -transparent");
        elseif strcmpi(grbType,"sgrb")
            export_fig(synSam.output.path + "/" + "sgrbDNDTObservedMerged.png", "-m4 -transparent");
        elseif strcmpi(grbType,"all")
            export_fig(synSam.output.path + "/" + "allDNDTObservedMerged.png", "-m4 -transparent");
        end
    end
hold off;