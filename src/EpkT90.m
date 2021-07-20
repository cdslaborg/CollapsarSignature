clear all;
close all;
filePath = mfilename('fullpath');
[currentDir,fileName,fileExt] = fileparts(filePath); cd(currentDir);
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath(genpath("../../../libmatlab"), "-begin");

% Figure Parameters
fontSize = 14;
lineWidth = 1.5;
figureColor = "white";
defaultColor = [0, 0.4470, 0.7410];
secondaryColor = '#FFA500'; %orange

% Simulation Options
freshRunEnabled = false; % this must be set to 'true' for first ever simulation. Thereafter, it can be set to 'false' to save time.
mergeBins = false; % merge adjacent bins with <5 events
saveNewImages = false; % export figure on or off. Must have export_fig installed.
skip = 28; % reduce synthetic lgrb dataset of 200,000 to ~1400, equivalent to BATSE; 'freshRunEnabled' must be 'true'
nbins = 50; % number of bins to bin synthetic data in; 'freshRunEnabled' must be 'true'
nruns = 10000; % number of detector simulation runs; 'freshRunEnabled' must be 'true'
matFileName = "EpkT90.mat";

if freshRunEnabled
    
    lgrb = importdata("..\in\syntheticSampleB10.csv");
    
    icol = struct();
    icol.logEpk = 6; % column index of logEpk
    icol.logT90 = 8; % column index of logDur
    icol.detProb = 10; % column index of detection probability
    
    synSam = struct();
    synSam.whole.Epk = exp(lgrb.data(:,icol.logEpk));
    synSam.whole.T90 = exp(lgrb.data(:,icol.logT90));
    synSam.whole.detProb = lgrb.data(:,icol.detProb);
    
    % Create bins
    synSam.whole.binEdges = getBinEdges(synSam.whole.T90, nbins);
    synSam.whole.binCenters = getBinCenters(synSam.whole.binEdges);
    
    % Determine BATSE-detectable events and reduce dataset size to be comparable to BATSE LGRB catalog. Repeat 'nruns' times.
    for i = 1:nruns
        % If detection probability is above a uniformly generated random number, observe event
        Mask = synSam.whole.detProb > unifrnd(0, 1, length(synSam.whole.detProb), 1);
        synSam.all.detectedVec(i).Epk = synSam.whole.Epk(Mask, 1);
        synSam.all.detectedVec(i).T90 = synSam.whole.T90(Mask, 1);
        synSam.all.detectedVec(i).Epk = synSam.all.detectedVec(i).Epk(1:skip:end);
        synSam.all.detectedVec(i).T90 = synSam.all.detectedVec(i).T90(1:skip:end);
        synSam.all.detectedVec(i).countsVec = histcounts(synSam.all.detectedVec(i).T90, synSam.whole.binEdges);
        
        % Cut the data set at the mean
        synSam.cut50.detectedVec(i).thresh = exp(mean(log(synSam.all.detectedVec(i).Epk)));
        Mask = synSam.all.detectedVec(i).Epk < synSam.cut50.detectedVec(i).thresh;
        synSam.cut50.detectedVec(i).Epk = synSam.all.detectedVec(i).Epk(Mask, 1);
        synSam.cut50.detectedVec(i).T90 = synSam.all.detectedVec(i).T90(Mask, 1);
        synSam.cut50.detectedVec(i).countsVec = histcounts(synSam.cut50.detectedVec(i).T90, synSam.whole.binEdges);
        
        % Cut the data set one std above mean
        synSam.cut84.detectedVec(i).thresh = exp(mean(log(synSam.all.detectedVec(i).Epk)) + std(log(synSam.all.detectedVec(i).Epk)));
        Mask = synSam.all.detectedVec(i).Epk < synSam.cut84.detectedVec(i).thresh;
        synSam.cut84.detectedVec(i).Epk = synSam.all.detectedVec(i).Epk(Mask, 1);
        synSam.cut84.detectedVec(i).T90 = synSam.all.detectedVec(i).T90(Mask, 1);
        synSam.cut84.detectedVec(i).countsVec = histcounts(synSam.cut84.detectedVec(i).T90, synSam.whole.binEdges);
        
        % Cut the data set one std below mean
        synSam.cut16.detectedVec(i).thresh = exp(mean(log(synSam.all.detectedVec(i).Epk)) - std(log(synSam.all.detectedVec(i).Epk)));
        Mask = synSam.all.detectedVec(i).Epk < synSam.cut16.detectedVec(i).thresh;
        synSam.cut16.detectedVec(i).Epk = synSam.all.detectedVec(i).Epk(Mask, 1);
        synSam.cut16.detectedVec(i).T90 = synSam.all.detectedVec(i).T90(Mask, 1);
        synSam.cut16.detectedVec(i).countsVec = histcounts(synSam.cut16.detectedVec(i).T90, synSam.whole.binEdges);
    end
    clear('lgrb', 'icol', 'Mask');
    
    % Grab a random set of detected data for Figure 1
    i = randi(nruns);
    synSam.all.detectedRand.Epk = synSam.all.detectedVec(i).Epk;
    synSam.all.detectedRand.T90 = synSam.all.detectedVec(i).T90;
    synSam.cut50.detectedRand.Epk = synSam.cut50.detectedVec(i).Epk;
    synSam.cut50.detectedRand.T90 = synSam.cut50.detectedVec(i).T90;
    synSam.cut84.detectedRand.Epk = synSam.cut84.detectedVec(i).Epk;
    synSam.cut84.detectedRand.T90 = synSam.cut84.detectedVec(i).T90;
    synSam.cut16.detectedRand.Epk = synSam.cut16.detectedVec(i).Epk;
    synSam.cut16.detectedRand.T90 = synSam.cut16.detectedVec(i).T90;
    synSam.cut50.detectedRand.thresh = synSam.cut50.detectedVec(i).thresh;
    synSam.cut84.detectedRand.thresh = synSam.cut84.detectedVec(i).thresh;
    synSam.cut16.detectedRand.thresh = synSam.cut16.detectedVec(i).thresh;
    
    for i = 1:nbins
        % Switch the positions of countsVec and detectedVec in the struct
        for j = 1:nruns
            synSam.all.countsVec(i).detectedVec(j) = synSam.all.detectedVec(j).countsVec(i);
            synSam.cut50.countsVec(i).detectedVec(j) = synSam.cut50.detectedVec(j).countsVec(i);
            synSam.cut84.countsVec(i).detectedVec(j) = synSam.cut84.detectedVec(j).countsVec(i);
            synSam.cut16.countsVec(i).detectedVec(j) = synSam.cut16.detectedVec(j).countsVec(i);
        end
        % Calculate the mean, the 95th percentile, and the 5th percentile
        synSam.all.percentile50.countsVec(i) = round(prctile(synSam.all.countsVec(i).detectedVec, 50));
        synSam.cut50.percentile50.countsVec(i) = round(prctile(synSam.cut50.countsVec(i).detectedVec, 50));
        synSam.cut84.percentile50.countsVec(i) = round(prctile(synSam.cut84.countsVec(i).detectedVec, 50));
        synSam.cut16.percentile50.countsVec(i) = round(prctile(synSam.cut16.countsVec(i).detectedVec, 50));
        %synSam.all.percentile95.countsVec(i) = round(prctile(synSam.all.countsVec(i).detectedVec, 95));
        %synSam.cut50.percentile95.countsVec(i) = round(prctile(synSam.cut50.countsVec(i).detectedVec, 95));
        %synSam.cut84.percentile95.countsVec(i) = round(prctile(synSam.cut84.countsVec(i).detectedVec, 95));
        %synSam.cut16.percentile95.countsVec(i) = round(prctile(synSam.cut16.countsVec(i).detectedVec, 95));
        %synSam.all.percentile05.countsVec(i) = round(prctile(synSam.all.countsVec(i).detectedVec, 5));
        %synSam.cut50.percentile05.countsVec(i) = round(prctile(synSam.cut50.countsVec(i).detectedVec, 5));
        %synSam.cut84.percentile05.countsVec(i) = round(prctile(synSam.cut84.countsVec(i).detectedVec, 5));
        %synSam.cut16.percentile05.countsVec(i) = round(prctile(synSam.cut16.countsVec(i).detectedVec, 5));
    end
    clear('i', 'j');
    synSam.all = rmfield(synSam.all, {'countsVec', 'detectedVec'});
    synSam.cut50 = rmfield(synSam.cut50, {'countsVec', 'detectedVec'});
    synSam.cut84 = rmfield(synSam.cut84, {'countsVec', 'detectedVec'});
    synSam.cut16 = rmfield(synSam.cut16, {'countsVec', 'detectedVec'});
    
    % Calculate ndata
    synSam.all.percentile50.ndata = sum(synSam.all.percentile50.countsVec);
    synSam.cut50.percentile50.ndata = sum(synSam.cut50.percentile50.countsVec);
    synSam.cut84.percentile50.ndata = sum(synSam.cut84.percentile50.countsVec);
    synSam.cut16.percentile50.ndata = sum(synSam.cut16.percentile50.countsVec);
    
    % Convert the dependent variable dN/dlog(T) to dN/dT by dividing it by T90
    synSam.all.percentile50.dNdT = synSam.all.percentile50.countsVec ./ synSam.whole.binCenters;
    synSam.cut50.percentile50.dNdT = synSam.cut50.percentile50.countsVec ./ synSam.whole.binCenters;
    synSam.cut84.percentile50.dNdT = synSam.cut84.percentile50.countsVec ./ synSam.whole.binCenters;
    synSam.cut16.percentile50.dNdT = synSam.cut16.percentile50.countsVec ./ synSam.whole.binCenters;
    
    % Save mat file
    synSam.output.path = "../out"; if ~isfolder(synSam.output.path); mkdir(synSam.output.path); end
    save(synSam.output.path + "/" + matFileName,"synSam");
else
    % Load mat file
    synSam.output.path = "../out";
    load(synSam.output.path + "/" + matFileName); % loads synSam object   
end

% Figure 1: Epk vs T90 of dataset
figure("color", figureColor); hold on; box on;
    scatter(synSam.all.detectedRand.T90, synSam.all.detectedRand.Epk, 5, defaultColor, 'filled', 'MarkerFaceAlpha', 0.25);
    scatter(synSam.cut84.detectedRand.T90, synSam.cut84.detectedRand.Epk, 5, defaultColor, 'filled', 'MarkerFaceAlpha', 0.25);
    scatter(synSam.cut50.detectedRand.T90, synSam.cut50.detectedRand.Epk, 5, defaultColor, 'filled', 'MarkerFaceAlpha', 0.25);
    scatter(synSam.cut16.detectedRand.T90, synSam.cut16.detectedRand.Epk, 5, defaultColor, 'filled');
    yline([synSam.cut84.detectedRand.thresh, synSam.cut50.detectedRand.thresh, synSam.cut16.detectedRand.thresh]);
    [ylower, yupper] = plotLimits(synSam.all.detectedRand.T90, 'log');
    xlim([ylower, yupper]);
    [ylower, yupper] = plotLimits(synSam.all.detectedRand.Epk, 'log');
    ylim([ylower, yupper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("E_{peak} [keV]", "interpreter", "tex", "fontsize", fontSize);
    title("Detectable LGRBs (n = " + synSam.all.percentile50.ndata + ")", "fontSize", fontSize+1);
hold off;

% Figure 2: Plot dN/dT of full detected dataset
figure("color", figureColor); hold on; box on;
    hist2stairs(synSam.all.percentile50.dNdT, synSam.whole.binEdges, defaultColor, lineWidth);
    if mergeBins
        newBins = binMergeStairPlot(synSam.all.percentile50.countsVec, synSam.whole.binEdges, secondaryColor, lineWidth);
        legend([synSam.all.percentile50.ndata + " detected LGRBs" + newline + "spanning " + nbins + " bins" ...
               , "Adjacent bins with <5 events merged," + newline + "resulting in " + newBins + " bins"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Full Detected Data (mean of " + nruns + " runs)", "fontSize", fontSize+1);
    else
        legend([synSam.all.percentile50.ndata + " detected LGRBs"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Full Detected Data (mean of " + nruns + " runs)", "fontSize", fontSize+1);
    end
    [xlower, xupper] = plotLimits(synSam.whole.binCenters, 'log', synSam.all.percentile50.dNdT);
    xlim([xlower, xupper]);
    [ylower, yupper] = plotLimits(synSam.all.percentile50.dNdT, 'log');
    ylim([ylower, yupper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
hold off;

% Figure 3: Plot dN/dT of dataset cut at mean
figure("color", figureColor); hold on; box on;
    hist2stairs(synSam.cut50.percentile50.dNdT, synSam.whole.binEdges, defaultColor, lineWidth);
    if mergeBins
        newBins = binMergeStairPlot(synSam.cut50.percentile50.countsVec, synSam.whole.binEdges, secondaryColor, lineWidth);
        legend([synSam.cut50.percentile50.ndata + " detected LGRBs" + newline + "spanning " + nbins + " bins" ...
               , "Adjacent bins with <5 events merged," + newline + "resulting in " + newBins + " bins"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean (mean of " + nruns + " runs)", "fontSize", fontSize+1);
    else
        legend([synSam.cut50.percentile50.ndata + " detected LGRBs"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean (mean of " + nruns + " runs)", "fontSize", fontSize+1);
    end
    xlim([xlower, xupper]);
    [ylower, yupper] = plotLimits(synSam.cut50.percentile50.dNdT, 'log');
    ylim([ylower, yupper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
hold off;

% Figure 4: Plot dN/dT of dataset cut at mean + std
figure("color", figureColor); hold on; box on;
    hist2stairs(synSam.cut84.percentile50.dNdT, synSam.whole.binEdges, defaultColor, lineWidth);
    if mergeBins
        newBins = binMergeStairPlot(synSam.cut84.percentile50.countsVec, synSam.whole.binEdges, secondaryColor, lineWidth);
        legend([synSam.cut84.percentile50.ndata + " detected LGRBs" + newline + "spanning " + nbins + " bins" ...
               , "Adjacent bins with <5 events merged," + newline + "resulting in " + newBins + " bins"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean + Std (mean of " + nruns + " runs)", "fontSize", fontSize+1);
    else
        legend([synSam.cut84.percentile50.ndata + " detected LGRBs"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean + Std (mean of " + nruns + " runs)", "fontSize", fontSize+1);
    end
    xlim([xlower, xupper]);
    [ylower, yupper] = plotLimits(synSam.cut84.percentile50.dNdT, 'log');
    ylim([ylower, yupper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
hold off;

% Figure 5: Plot dN/dT of dataset cut at mean - std
figure("color", figureColor); hold on; box on;
    hist2stairs(synSam.cut16.percentile50.dNdT, synSam.whole.binEdges, defaultColor, lineWidth);
    if mergeBins
        newBins = binMergeStairPlot(synSam.cut16.percentile50.countsVec, synSam.whole.binEdges, secondaryColor, lineWidth);
        legend([synSam.cut16.percentile50.ndata + " detected LGRBs" + newline + "spanning " + nbins + " bins" ...
               , "Adjacent bins with <5 events merged," + newline + "resulting in " + newBins + " bins"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean - Std (mean of " + nruns + " runs)", "fontSize", fontSize+1);
    else
        legend([synSam.cut16.percentile50.ndata + " detected LGRBs"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean - Std (mean of " + nruns + " runs)", "fontSize", fontSize+1);
    end
    xlim([xlower, xupper]);
    [ylower, yupper] = plotLimits(synSam.cut16.percentile50.dNdT, 'log');
    ylim([ylower, yupper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
hold off;