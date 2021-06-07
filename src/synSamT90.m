clear all;
close all;
filePath = mfilename('fullpath');
[currentDir,fileName,fileExt] = fileparts(filePath); cd(currentDir);
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath(genpath("../../../libmatlab"),"-begin");

% Figure Parameters
fontSize = 14;
lineWidth = 1.5;
figureColor = "white";
secondaryColor = "magenta";

% Simulation Options
grbType = "all"; % choose between "lgrb", "sgrb", and "all"
freshRunEnabled = true; % this must be set to 'true' for first ever simulation. Thereafter, it can be set to 'false' to save time.
skip.lgrb = 28; % reduce synthetic lgrb dataset of 200,000 to ~1400, equivalent to BATSE; 'freshRunEnabled' must be 'true'
skip.sgrb = 10; % reduce synthetic sgrb dataset of 18,848 to ~600, equivalent to BATSE; 'freshRunEnabled' must be 'true'
numBins = 40; % number of bins to bin data in; 'freshRunEnabled' must be 'true'
numRuns = 10000; % number of detector simulation runs; 'freshRunEnabled' must be 'true'
matFileName = "synSamT90.mat";

if ~any(strcmpi(grbType,["lgrb","sgrb","all"])); error("'grbType' is not an allowable value"); end
if freshRunEnabled
    
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
    synSam.lgrb.binEdges = exp(linspace(log(min(synSam.lgrb.T90)),log(max(synSam.lgrb.T90)),numBins+1));
    synSam.lgrb.binCenters = zeros(1,numBins);
    synSam.sgrb.binEdges = exp(linspace(log(min(synSam.sgrb.T90)),log(max(synSam.sgrb.T90)),numBins+1));
    synSam.sgrb.binCenters = zeros(1,numBins);
    synSam.all.binEdges = exp(linspace(log(min(synSam.sgrb.T90)),log(max(synSam.lgrb.T90)),numBins+1));
    synSam.all.binCenters = zeros(1,numBins);
    for i = 1:numBins
        synSam.lgrb.binCenters(i) = 0.5*(synSam.lgrb.binEdges(i) + synSam.lgrb.binEdges(i+1));
        synSam.sgrb.binCenters(i) = 0.5*(synSam.sgrb.binEdges(i) + synSam.sgrb.binEdges(i+1));
        synSam.all.binCenters(i) = 0.5*(synSam.all.binEdges(i) + synSam.all.binEdges(i+1));
    end
    synSam.lgrb.whole.ndata = length(synSam.lgrb.T90);
    synSam.lgrb.whole.countsVec = histcounts(synSam.lgrb.T90,synSam.lgrb.binEdges);
    synSam.sgrb.whole.ndata = length(synSam.sgrb.T90);
    synSam.sgrb.whole.countsVec = histcounts(synSam.sgrb.T90,synSam.sgrb.binEdges);
    synSam.all.whole.ndata = length([synSam.lgrb.T90;synSam.sgrb.T90]);
    synSam.all.whole.countsVec = histcounts([synSam.lgrb.T90;synSam.sgrb.T90],synSam.all.binEdges);
    % Normalize the ratio of LGRBs/SGRBs to the observed ratio for BATSE in the combined, pre-detected dataset. Does not affect any other computations.
    synSam.all.whole.normalized.ndata = length([synSam.lgrb.T90(1:round(length(synSam.sgrb.T90)*1366/600));synSam.sgrb.T90]);
    synSam.all.whole.normalized.countsVec = histcounts([synSam.lgrb.T90(1:round(length(synSam.sgrb.T90)*1366/600));synSam.sgrb.T90],synSam.all.binEdges);
    for i = 1:numRuns
        Mask.lgrb = lgrb.data(:,icol.lgrb.detProb) > unifrnd(0,1,length(lgrb.data),1);
        Mask.sgrb = sgrb.data(:,icol.sgrb.detProb) > unifrnd(0,1,length(sgrb.data),1);
        T90.lgrb = exp(lgrb.data(Mask.lgrb,icol.lgrb.logT90));  % choose those events that would be detected
        T90.sgrb = sgrb.data(Mask.sgrb,icol.sgrb.T90);          % choose those events that would be detected
        T90.lgrb = T90.lgrb(1:skip.lgrb:end);   % reduce the detected dataset
        T90.sgrb = T90.sgrb(1:skip.sgrb:end);   % reduce the detected dataset
        synSam.lgrb.reducedVec(i).ndata = length(T90.lgrb);
        synSam.lgrb.reducedVec(i).countsVec = histcounts(T90.lgrb,synSam.lgrb.binEdges);
        synSam.sgrb.reducedVec(i).ndata = length(T90.sgrb);
        synSam.sgrb.reducedVec(i).countsVec = histcounts(T90.sgrb,synSam.sgrb.binEdges);
        synSam.all.reducedVec(i).ndata = length([T90.lgrb;T90.sgrb]);
        synSam.all.reducedVec(i).countsVec = histcounts([T90.lgrb;T90.sgrb],synSam.all.binEdges);
    end
    clear('lgrb','sgrb','icol','Mask');
    
    for i = 1:numBins
        for j = 1:numRuns
            synSam.lgrb.countsVec(i).reducedVec(j) = synSam.lgrb.reducedVec(j).countsVec(i);
            synSam.sgrb.countsVec(i).reducedVec(j) = synSam.sgrb.reducedVec(j).countsVec(i);
            synSam.all.countsVec(i).reducedVec(j) = synSam.all.reducedVec(j).countsVec(i);
        end
        synSam.lgrb.reduced.percentile95.countsVec(i) = prctile(synSam.lgrb.countsVec(i).reducedVec,95);
        synSam.lgrb.reduced.percentile50.countsVec(i) = prctile(synSam.lgrb.countsVec(i).reducedVec,50);
        synSam.lgrb.reduced.percentile05.countsVec(i) = prctile(synSam.lgrb.countsVec(i).reducedVec,5);
        synSam.sgrb.reduced.percentile95.countsVec(i) = prctile(synSam.sgrb.countsVec(i).reducedVec,95);
        synSam.sgrb.reduced.percentile50.countsVec(i) = prctile(synSam.sgrb.countsVec(i).reducedVec,50);
        synSam.sgrb.reduced.percentile05.countsVec(i) = prctile(synSam.sgrb.countsVec(i).reducedVec,5);
        synSam.all.reduced.percentile95.countsVec(i) = prctile(synSam.all.countsVec(i).reducedVec,95);
        synSam.all.reduced.percentile50.countsVec(i) = prctile(synSam.all.countsVec(i).reducedVec,50);
        synSam.all.reduced.percentile05.countsVec(i) = prctile(synSam.all.countsVec(i).reducedVec,5);
    end
    clear('i','j');
    synSam.lgrb.reduced.mean.ndata = round(mean([synSam.lgrb.reducedVec(:).ndata]));
    synSam.sgrb.reduced.mean.ndata = round(mean([synSam.sgrb.reducedVec(:).ndata]));
    synSam.all.reduced.mean.ndata = round(mean([synSam.all.reducedVec(:).ndata]));
    
    synSam.output.path = "../out"; if ~isfolder(synSam.output.path); mkdir(synSam.output.path); end
    save(synSam.output.path + "/" + matFileName,"synSam");

else
    
    synSam.output.path = "../out";
    load(synSam.output.path + "/" + matFileName); % loads synSam object
    
end

% Figure 1: Histogram for whole dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        histogram('BinEdges', synSam.lgrb.binEdges, 'BinCounts', synSam.lgrb.whole.countsVec);
        title("Histogram for whole LGRB dataset n = " + synSam.lgrb.whole.ndata, "fontsize", fontSize);
    elseif strcmpi(grbType,"sgrb")
        histogram('BinEdges', synSam.sgrb.binEdges, 'BinCounts', synSam.sgrb.whole.countsVec);
        title("Histogram for whole SGRB dataset n = " + synSam.sgrb.whole.ndata, "fontsize", fontSize);
        ylim([1,10^4]);
    elseif strcmpi(grbType,"all")
        histogram('BinEdges', synSam.all.binEdges, 'BinCounts', synSam.all.whole.countsVec);
        histogram('BinEdges', synSam.all.binEdges, 'BinCounts', synSam.all.whole.normalized.countsVec, 'FaceColor', secondaryColor);
        %title("Histogram for whole GRB dataset n = " + synSam.all.whole.ndata, "fontsize", fontSize);
        title("Histogram for whole GRB dataset", "fontsize", fontSize);
        legend(["Whole dataset, n = " + synSam.all.whole.ndata, "LGRB/SGRB ratio normalized to BATSE, n = " + synSam.all.whole.normalized.ndata] ...
                , "interpreter", "tex", "location", "northwest");
        ylim([1,10^6]);
    end
    xlabel("T_{90}", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN/dlog(T_{90})", "interpreter", "tex", "fontsize", fontSize);
    set(gca, 'xscale', 'log', 'yscale', 'log');
hold off;

lgrb.dNdT = synSam.lgrb.whole.countsVec ./ synSam.lgrb.binCenters;
sgrb.dNdT = synSam.sgrb.whole.countsVec ./ synSam.sgrb.binCenters;
all.dNdT = synSam.all.whole.countsVec ./ synSam.all.binCenters;
all.dNdTnorm = synSam.all.whole.normalized.countsVec ./ synSam.all.binCenters;

% Figure 2: Plot dN/dT of whole dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        stairs(synSam.lgrb.binEdges(1:end-1), lgrb.dNdT, 'lineWidth', lineWidth);
        title("Adjusted LGRB dataset of n = " + synSam.lgrb.whole.ndata + " with original binning", "fontsize", fontSize);
    elseif strcmpi(grbType,"sgrb")
        stairs(synSam.sgrb.binEdges(1:end-1), sgrb.dNdT, 'lineWidth', lineWidth);
        title("Adjusted SGRB dataset of n = " + synSam.sgrb.whole.ndata + " with original binning", "fontsize", fontSize);
        ylim([10^-2,10^5]);
    elseif strcmpi(grbType,"all")
        stairs(synSam.all.binEdges(1:end-1), all.dNdT, 'lineWidth', lineWidth);
        stairs(synSam.all.binEdges(1:end-1), all.dNdTnorm, 'lineWidth', lineWidth, 'color', secondaryColor);
        %title("Adjusted GRB dataset of n = " + synSam.all.whole.ndata + " with original binning", "fontsize", fontSize);
        title("Adjusted GRB dataset with original binning", "fontsize", fontSize);
        legend(["Whole dataset, n = " + synSam.all.whole.ndata, "LGRB/SGRB ratio normalized to BATSE, n = " + synSam.all.whole.normalized.ndata] ...
                , "interpreter", "tex", "location", "southwest");
        ylim([10^-4,10^5]);
    end
    xlabel("T_{90}", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN/dT_{90}", "interpreter", "tex", "fontsize", fontSize);
    set(gca, 'xscale', 'log', 'yscale', 'log');
hold off;

lgrb.dNdT95th = synSam.lgrb.reduced.percentile95.countsVec ./ synSam.lgrb.binCenters;
lgrb.dNdT05th = synSam.lgrb.reduced.percentile05.countsVec ./ synSam.lgrb.binCenters;
sgrb.dNdT95th = synSam.sgrb.reduced.percentile95.countsVec ./ synSam.sgrb.binCenters;
sgrb.dNdT05th = synSam.sgrb.reduced.percentile05.countsVec ./ synSam.sgrb.binCenters;
all.dNdT50th = synSam.all.reduced.percentile50.countsVec ./ synSam.all.binCenters;
all.dNdT95th = synSam.all.reduced.percentile95.countsVec ./ synSam.all.binCenters;
all.dNdT05th = synSam.all.reduced.percentile05.countsVec ./ synSam.all.binCenters;

lgrb.zeroIndicies = find(~synSam.lgrb.reduced.percentile50.countsVec);
synSam.lgrb.reduced.percentile50.countsVecTrimmed = synSam.lgrb.reduced.percentile50.countsVec;
synSam.lgrb.reduced.percentile50.countsVecTrimmed(lgrb.zeroIndicies) = [];
synSam.lgrb.binCentersTrimmed = synSam.lgrb.binCenters;
synSam.lgrb.binCentersTrimmed(lgrb.zeroIndicies) = [];
synSam.lgrb.binEdgesTrimmed = synSam.lgrb.binEdges;
synSam.lgrb.binEdgesTrimmed(lgrb.zeroIndicies) = [];
lgrb.fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[130,3,1]);
lgrb.ft = fittype('a*exp(-((log(x)-b)/c)^2)', 'options', lgrb.fo);
lgrb.f = fit(synSam.lgrb.binCentersTrimmed.', synSam.lgrb.reduced.percentile50.countsVecTrimmed.', lgrb.ft);
lgrb.dNdT50th = synSam.lgrb.reduced.percentile50.countsVecTrimmed ./ synSam.lgrb.binCentersTrimmed;
lgrb.fdNdT = lgrb.f(synSam.lgrb.binCentersTrimmed).' ./ synSam.lgrb.binCentersTrimmed;

sgrb.zeroIndicies = find(~synSam.sgrb.reduced.percentile50.countsVec);
synSam.sgrb.reduced.percentile50.countsVecTrimmed = synSam.sgrb.reduced.percentile50.countsVec;
synSam.sgrb.reduced.percentile50.countsVecTrimmed(sgrb.zeroIndicies) = [];
synSam.sgrb.binCentersTrimmed = synSam.sgrb.binCenters;
synSam.sgrb.binCentersTrimmed(sgrb.zeroIndicies) = [];
synSam.sgrb.binEdgesTrimmed = synSam.sgrb.binEdges;
synSam.sgrb.binEdgesTrimmed(sgrb.zeroIndicies) = [];
sgrb.fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,1,1]);
sgrb.ft = fittype('a*exp(-((log(x)-b)/c)^2)', 'options', sgrb.fo);
sgrb.f = fit(synSam.sgrb.binCentersTrimmed.', synSam.sgrb.reduced.percentile50.countsVecTrimmed.', sgrb.ft);
sgrb.dNdT50th = synSam.sgrb.reduced.percentile50.countsVecTrimmed ./ synSam.sgrb.binCentersTrimmed;
sgrb.fdNdT = sgrb.f(synSam.sgrb.binCentersTrimmed).' ./ synSam.sgrb.binCentersTrimmed;

% Figure 3: Plot dN/dT of reduced dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        stairs(synSam.lgrb.binEdges(1:end-1), lgrb.dNdT95th, 'color', 'red');
        stairs(synSam.lgrb.binEdges(1:end-1), lgrb.dNdT05th, 'color', 'red');
        stairs(synSam.lgrb.binEdges(1:end-1), lgrb.dNdT50th, 'color', [0, 0.4470, 0.7410], 'lineWidth', lineWidth);
        title("Adjusted LGRB dataset of n = " + synSam.lgrb.reduced.mean.ndata + " with original binning", "fontsize", fontSize);
        legend(["", "90% Confidence Interval", "Mean of " + numRuns + " runs"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize);
    elseif strcmpi(grbType,"sgrb")
        stairs(synSam.sgrb.binEdges(1:end-1), sgrb.dNdT95th, 'color', 'red');
        stairs(synSam.sgrb.binEdges(1:end-1), sgrb.dNdT05th, 'color', 'red');
        stairs(synSam.sgrb.binEdges(1:end-1), sgrb.dNdT50th, 'color', [0, 0.4470, 0.7410], 'lineWidth', lineWidth);
        title("Adjusted SGRB dataset of n = " + synSam.sgrb.reduced.mean.ndata + " with original binning", "fontsize", fontSize);
        legend(["", "90% Confidence Interval", "Mean of " + numRuns + " runs"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize);
    elseif strcmpi(grbType,"all")
        %stairs(synSam.all.binEdges(1:end-1), all.dNdT95th, 'color', 'red');
        %stairs(synSam.all.binEdges(1:end-1), all.dNdT05th, 'color', 'red');
        %stairs(synSam.all.binEdges(1:end-1), all.dNdT50th, 'color', [0.5, 0.5, 0.5], 'lineWidth', lineWidth-0.5);
        plot(synSam.sgrb.binCentersTrimmed, sgrb.fdNdT, synSam.lgrb.binCentersTrimmed, lgrb.fdNdT, 'color', [0.5, 0.5, 0.5], 'lineWidth', lineWidth);
        plot( cat(2, synSam.sgrb.binCentersTrimmed, synSam.lgrb.binCentersTrimmed) ...
            , cat(2, sgrb.fdNdT, lgrb.fdNdT), 'color', [0.5, 0.5, 0.5], 'lineWidth', lineWidth);
        plot(synSam.lgrb.binCenters, lgrb.dNdT95th, ':', 'color', 'green', 'lineWidth', lineWidth);
        plot(synSam.lgrb.binCenters, lgrb.dNdT05th, ':', 'color', 'green', 'lineWidth', lineWidth);
        stairs(synSam.lgrb.binEdgesTrimmed(1:end-1), lgrb.dNdT50th, 'color', 'red', 'lineWidth', lineWidth-0.5);
        plot(synSam.lgrb.binCentersTrimmed, lgrb.fdNdT, 'color', 'red', 'lineWidth', lineWidth);
        plot(synSam.sgrb.binCenters, sgrb.dNdT95th, ':', 'color', 'green', 'lineWidth', lineWidth);
        plot(synSam.sgrb.binCenters, sgrb.dNdT05th, ':', 'color', 'green', 'lineWidth', lineWidth);
        stairs(synSam.sgrb.binEdgesTrimmed(1:end-1), sgrb.dNdT50th, 'color', 'blue', 'lineWidth', lineWidth-0.5);
        plot(synSam.sgrb.binCentersTrimmed, sgrb.fdNdT, 'color', 'blue', 'lineWidth', lineWidth);
        title("Adjusted GRB dataset of n = " + synSam.all.reduced.mean.ndata + " with original binning", "fontsize", fontSize);
        xlim([10^-2, 2*10^3]);
    end
    xlabel("T_{90}", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN/dT_{90}", "interpreter", "tex", "fontsize", fontSize);
    set(gca, 'xscale', 'log', 'yscale', 'log');
hold off;

% Figure 4: Combine adjacent bins with fewer than 5 events and plot dN/dT of reduced dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        binReductionStairPlot(synSam.lgrb.reduced.percentile95.countsVec, synSam.lgrb.binEdges, 'red');
        binReductionStairPlot(synSam.lgrb.reduced.percentile05.countsVec, synSam.lgrb.binEdges, 'red');
        nbins = binReductionStairPlot(synSam.lgrb.reduced.percentile50.countsVec, synSam.lgrb.binEdges, [0, 0.4470, 0.7410], lineWidth);
        title("Adjusted LGRB dataset of n = " + synSam.lgrb.reduced.mean.ndata + " in " + nbins + " bins", "fontsize", fontSize);
    elseif strcmpi(grbType,"sgrb")
        binReductionStairPlot(synSam.sgrb.reduced.percentile95.countsVec, synSam.sgrb.binEdges, 'red');
        binReductionStairPlot(synSam.sgrb.reduced.percentile05.countsVec, synSam.sgrb.binEdges, 'red');
        nbins = binReductionStairPlot(synSam.sgrb.reduced.percentile50.countsVec, synSam.sgrb.binEdges, [0, 0.4470, 0.7410], lineWidth);
        title("Adjusted SGRB dataset of n = " + synSam.sgrb.reduced.mean.ndata + " in " + nbins + " bins", "fontsize", fontSize);
    elseif strcmpi(grbType,"all")
        binReductionStairPlot(synSam.all.reduced.percentile95.countsVec, synSam.all.binEdges, 'red');
        binReductionStairPlot(synSam.all.reduced.percentile05.countsVec, synSam.all.binEdges, 'red');
        nbins = binReductionStairPlot(synSam.all.reduced.percentile50.countsVec, synSam.all.binEdges, [0, 0.4470, 0.7410], lineWidth);
        title("Adjusted GRB dataset of n = " + synSam.all.reduced.mean.ndata + " in " + nbins + " bins", "fontsize", fontSize);
    end
    xlabel("T_{90}", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN/dT_{90}", "interpreter", "tex", "fontsize", fontSize);
    legend(["", "90% Confidence Interval", "Mean of " + numRuns + " runs"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize);
    set(gca, 'xscale', 'log', 'yscale', 'log');
hold off;

% Figure 5: Histogram for reduced dataset of all
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', synSam.all.binEdges, 'BinCounts', synSam.all.reduced.percentile50.countsVec);
    xline(3, 'color', 'green');
    title("Histogram for reduced GRB dataset", "fontsize", fontSize);
    xlim([10^-2,10^3]);
    xlabel("T_{90}", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN/dlog(T_{90})", "interpreter", "tex", "fontsize", fontSize);
    set(gca, 'xscale', 'log');
hold off;
