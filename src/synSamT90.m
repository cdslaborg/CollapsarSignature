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
defaultColor = [0, 0.4470, 0.7410];
secondaryColor = '#FFA500'; %orange
sgrbColor = "blue";
lgrbColor = "red";
allColor = "magenta";
confidenceColor = [0.5, 0.5, 0.5];

% Simulation Options
grbType = "all"; % choose between "lgrb", "sgrb", and "all"
freshRunEnabled = false; % this must be set to 'true' for first ever simulation. Thereafter, it can be set to 'false' to save time.
skip.lgrb = 28; % reduce synthetic lgrb dataset of 200,000 to ~1400, equivalent to BATSE; 'freshRunEnabled' must be 'true'
skip.sgrb = 10; % reduce synthetic sgrb dataset of 18,848 to ~600, equivalent to BATSE; 'freshRunEnabled' must be 'true'
numBins = 52; % number of bins to bin data in; 'freshRunEnabled' must be 'true'
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
    
    % Trim zeros off 
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
    
    synSam.output.path = "../out"; if ~isfolder(synSam.output.path); mkdir(synSam.output.path); end
    save(synSam.output.path + "/" + matFileName,"synSam");

else
    
    synSam.output.path = "../out";
    load(synSam.output.path + "/" + matFileName); % loads synSam object
    
end

% Figure 1: Histogram for whole dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        histogram('BinEdges', synSam.lgrb.whole.binEdges, 'BinCounts', synSam.lgrb.whole.countsVec);
        title("Histogram for whole LGRB dataset n = " + synSam.lgrb.whole.ndata + " in " + length(synSam.lgrb.whole.countsVec) + " bins");
    elseif strcmpi(grbType,"sgrb")
        histogram('BinEdges',synSam.sgrb.whole.binEdges, 'BinCounts', synSam.sgrb.whole.countsVec);
        title("Histogram for whole SGRB dataset n = " + synSam.sgrb.whole.ndata + " in " + length(synSam.sgrb.whole.countsVec) + " bins");
    elseif strcmpi(grbType,"all")
        histogram('BinEdges', synSam.all.binEdges, 'BinCounts', synSam.all.whole.countsVec);
        histogram('BinEdges', synSam.all.whole.normalized.binEdges, 'BinCounts', synSam.all.whole.normalized.countsVec, 'FaceColor', secondaryColor);
        %title("Histogram for whole GRB dataset n = " + synSam.all.whole.ndata, "fontsize", fontSize);
        title("Histogram for whole GRB dataset in " + length(synSam.all.whole.countsVec) + " bins");
        legend(["Whole dataset, n = " + synSam.all.whole.ndata, "LGRB/SGRB ratio normalized to BATSE, n = " + synSam.all.whole.normalized.ndata] ...
                , "interpreter", "tex", "location", "northwest");
        ylim([1,10^6]);
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dlog(T_{90})", "interpreter", "tex", "fontsize", fontSize);
hold off;

lgrb.dNdT = synSam.lgrb.whole.countsVec ./ synSam.lgrb.whole.binCenters;
sgrb.dNdT = synSam.sgrb.whole.countsVec ./ synSam.sgrb.whole.binCenters;
all.dNdT = synSam.all.whole.countsVec ./ synSam.all.binCenters;
all.dNdTnorm = synSam.all.whole.normalized.countsVec ./ synSam.all.whole.normalized.binCenters;

% Figure 2: Plot dN/dT of whole dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        hist2stairs(lgrb.dNdT, synSam.lgrb.whole.binEdges, defaultColor, lineWidth);
        title("Adjusted LGRB dataset of n = " + synSam.lgrb.whole.ndata + " in " + length(lgrb.dNdT) + " bins");
    elseif strcmpi(grbType,"sgrb")
        hist2stairs(sgrb.dNdT, synSam.sgrb.whole.binEdges, defaultColor, lineWidth);
        title("Adjusted SGRB dataset of n = " + synSam.sgrb.whole.ndata + " in " + length(sgrb.dNdT) + " bins");
        ylim([10^-2, 10^5]);
    elseif strcmpi(grbType,"all")
        hist2stairs(all.dNdT, synSam.all.binEdges, defaultColor, lineWidth);
        hist2stairs(all.dNdTnorm, synSam.all.whole.normalized.binEdges, secondaryColor, lineWidth);
        title("Adjusted GRB dataset with original binning");
        legend(["Whole dataset, n = " + synSam.all.whole.ndata, "LGRB/SGRB ratio normalized to BATSE, n = " + synSam.all.whole.normalized.ndata] ...
                , "interpreter", "tex", "location", "southwest");
        ylim([10^-4, 10^5]);
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
hold off;

lgrb.dNdT95th = synSam.lgrb.reduced.percentile95.countsVec ./ synSam.lgrb.reduced.percentile95.binCenters;
lgrb.dNdT05th = synSam.lgrb.reduced.percentile05.countsVec ./ synSam.lgrb.reduced.percentile05.binCenters;
sgrb.dNdT95th = synSam.sgrb.reduced.percentile95.countsVec ./ synSam.sgrb.reduced.percentile95.binCenters;
sgrb.dNdT05th = synSam.sgrb.reduced.percentile05.countsVec ./ synSam.sgrb.reduced.percentile05.binCenters;

lgrb.fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[130,3,1]);
lgrb.ft = fittype('a*exp(-((log(x)-b)/c)^2)', 'options', lgrb.fo);
lgrb.f = fit(synSam.lgrb.reduced.percentile50.binCenters.', synSam.lgrb.reduced.percentile50.countsVec.', lgrb.ft);
lgrb.dNdT50th = synSam.lgrb.reduced.percentile50.countsVec ./ synSam.lgrb.reduced.percentile50.binCenters;
lgrb.fdNdT = lgrb.f(synSam.lgrb.reduced.percentile50.binCenters).' ./ synSam.lgrb.reduced.percentile50.binCenters;

sgrb.fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,1,1]);
sgrb.ft = fittype('a*exp(-((log(x)-b)/c)^2)', 'options', sgrb.fo);
sgrb.f = fit(synSam.sgrb.reduced.percentile50.binCenters.', synSam.sgrb.reduced.percentile50.countsVec.', sgrb.ft);
sgrb.dNdT50th = synSam.sgrb.reduced.percentile50.countsVec ./ synSam.sgrb.reduced.percentile50.binCenters;
sgrb.fdNdT = sgrb.f(synSam.sgrb.reduced.percentile50.binCenters).' ./ synSam.sgrb.reduced.percentile50.binCenters;

all.fdNdT = lgrb.f(synSam.all.reduced.percentile50.binCenters).' ./ synSam.all.reduced.percentile50.binCenters ...
          + sgrb.f(synSam.all.reduced.percentile50.binCenters).' ./ synSam.all.reduced.percentile50.binCenters;

%xq = exp(linspace(log(synSam.sgrb.binCentersTrimmed(1)),log(synSam.lgrb.binCentersTrimmed(end)),numBins+1));
%vq = interp1(synSam.sgrb.binCenters, sgrb.dNdT95th, xq, 'spline');

% Figure 3: Plot dN/dT of reduced dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        plot(synSam.lgrb.reduced.percentile95.binCenters, lgrb.dNdT95th, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        plot(synSam.lgrb.reduced.percentile05.binCenters, lgrb.dNdT05th, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        hist2stairs(lgrb.dNdT50th, synSam.lgrb.reduced.percentile50.binEdges, lgrbColor, lineWidth);
        title("Adjusted LGRB dataset of n = " + synSam.lgrb.reduced.mean.ndata + " in " + length(lgrb.dNdT50th) + " bins");
        legend(["", "90% confidence Interval", "Mean of " + numRuns + " runs"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize);
    elseif strcmpi(grbType,"sgrb")
        plot(synSam.sgrb.reduced.percentile95.binCenters, sgrb.dNdT95th, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        plot(synSam.sgrb.reduced.percentile05.binCenters, sgrb.dNdT05th, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        hist2stairs(sgrb.dNdT50th, synSam.sgrb.reduced.percentile50.binEdges, sgrbColor, lineWidth);
        title("Adjusted SGRB dataset of n = " + synSam.sgrb.reduced.mean.ndata + " in " + length(sgrb.dNdT50th) + " bins");
        legend(["", "90% confidence Interval", "Mean of " + numRuns + " runs"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize);
    elseif strcmpi(grbType,"all")
        plot(synSam.all.reduced.percentile50.binCenters, all.fdNdT, 'color', 'magenta', 'lineWidth', lineWidth+1);
        hist2stairs(lgrb.dNdT50th, synSam.lgrb.reduced.percentile50.binEdges, lgrbColor, lineWidth-0.5);
        plot(synSam.lgrb.reduced.percentile50.binCenters, lgrb.fdNdT, 'color', lgrbColor, 'lineWidth', lineWidth);
        hist2stairs(sgrb.dNdT50th, synSam.sgrb.reduced.percentile50.binEdges, sgrbColor, lineWidth-0.5);
        plot(synSam.sgrb.reduced.percentile50.binCenters, sgrb.fdNdT, 'color', sgrbColor, 'lineWidth', lineWidth);
        plot(synSam.lgrb.reduced.percentile95.binCenters, lgrb.dNdT95th, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        plot(synSam.lgrb.reduced.percentile05.binCenters, lgrb.dNdT05th, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        plot(synSam.sgrb.reduced.percentile95.binCenters, sgrb.dNdT95th, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        plot(synSam.sgrb.reduced.percentile05.binCenters, sgrb.dNdT05th, ':', 'color', confidenceColor, 'lineWidth', lineWidth);
        title("Adjusted GRB dataset of n = " + synSam.all.reduced.mean.ndata + " in " + length(all.fdNdT) + " bins");
        xlim([10^-2, 2*10^3]);
        legend([synSam.all.reduced.mean.ndata + " synthetic GRBs", "", synSam.lgrb.reduced.mean.ndata + " synthetic LGRBs" ...
               , "", synSam.sgrb.reduced.mean.ndata + " synthetic SGRBs", "90% confidence Interval"], "interpreter", "tex" ...
               , "location", "southwest", "fontSize", fontSize-2);
        legend('boxoff');
        ylim([10^-3,10^3]);
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
hold off;

% Figure 4: Combine adjacent bins with fewer than 5 events and plot dN/dT of reduced dataset
figure("color", figureColor); hold on; box on;
    if strcmpi(grbType,"lgrb")
        plot(synSam.lgrb.reduced.percentile50.binCenters, lgrb.fdNdT, 'color', lgrbColor, 'lineWidth', lineWidth+1);
        nbins = binMergeStairPlot(synSam.lgrb.reduced.percentile50.countsVec, synSam.lgrb.reduced.percentile50.binEdges, defaultColor, lineWidth);
        title("Adjusted LGRB dataset of n = " + synSam.lgrb.reduced.mean.ndata + " in " + nbins + " bins");
        legend(["BATSE " + synSam.lgrb.reduced.mean.ndata + " LGRB T_{90} dist.", "Adjacent bins with <5 events merged"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-2);
    elseif strcmpi(grbType,"sgrb")
        plot(synSam.sgrb.reduced.percentile50.binCenters, sgrb.fdNdT, 'color', sgrbColor, 'lineWidth', lineWidth+1);
        nbins = binMergeStairPlot(synSam.sgrb.reduced.percentile50.countsVec, synSam.sgrb.reduced.percentile50.binEdges, defaultColor, lineWidth);
        title("Adjusted SGRB dataset of n = " + synSam.sgrb.reduced.mean.ndata + " in " + nbins + " bins");
        legend(["BATSE " + synSam.sgrb.reduced.mean.ndata + " SGRB T_{90} dist.", "Adjacent bins with <5 events merged"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-2);
    elseif strcmpi(grbType,"all")
        plot(synSam.all.reduced.percentile50.binCenters, all.fdNdT, 'color', allColor, 'lineWidth', lineWidth+1);
        nbins = binMergeStairPlot(synSam.all.reduced.percentile50.countsVec, synSam.all.reduced.percentile50.binEdges, defaultColor, lineWidth);
        title("Adjusted GRB dataset of n = " + synSam.all.reduced.mean.ndata + " in " + nbins + " bins");
        legend([synSam.all.reduced.mean.ndata + " synthetic GRBs", "Adjacent bins with <5 events merged"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-2);
        xlim([10^-2,2*10^3]);
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
hold off;

% Figure 5: Histogram for reduced dataset of all
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', synSam.all.reduced.percentile50.binEdges, 'BinCounts', synSam.all.reduced.percentile50.countsVec);
    xline(3, 'color', 'green');
    title("Histogram for reduced GRB dataset", "fontsize", fontSize);
    xlim([10^-2,10^3]);
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dlog(T_{90})", "interpreter", "tex", "fontsize", fontSize);
hold off;
