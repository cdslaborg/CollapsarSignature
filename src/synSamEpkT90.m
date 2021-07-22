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

% Simulation Options
freshRunEnabled = false; % this must be set to 'true' for first ever simulation. Thereafter, it can be set to 'false' to save time.
mergeBins = false; % merge adjacent bins with <5 events
saveNewImages = false; % export figure on or off. Must have export_fig installed.
skip.lgrb = 28; % reduce synthetic lgrb dataset of 200,000 to ~1400, equivalent to BATSE; 'freshRunEnabled' must be 'true'
skip.sgrb = 10; % reduce synthetic sgrb dataset of 18,848 to ~600, equivalent to BATSE; 'freshRunEnabled' must be 'true'
nbins = 50; % number of bins to bin synthetic data in; 'freshRunEnabled' must be 'true'
nruns = 10000; % number of detector simulation runs; 'freshRunEnabled' must be 'true'
matFileName = "synSamEpkT90.mat";

if freshRunEnabled
    
    lgrb = importdata("..\in\syntheticSampleB10.csv");
    sgrb = importdata("..\in\SyntheticSample_All.txt");
    
    icol = struct();
    icol.lgrb.logEpk = 6; % column index of logEpk
    icol.lgrb.logT90 = 8; % column index of logDur
    icol.lgrb.detProb = 10; % column index of detection probability
    icol.sgrb.Epk = 7; % column index of Epk
    icol.sgrb.T90 = 9; % column index of T90
    icol.sgrb.detProb = 23; % column index of detection probability
    
    synSam = struct();
    synSam.lgrb.Epk = exp(lgrb.data(:,icol.lgrb.logEpk));
    synSam.lgrb.T90 = exp(lgrb.data(:,icol.lgrb.logT90));
    synSam.lgrb.detProb = lgrb.data(:,icol.lgrb.detProb);
    synSam.sgrb.Epk = sgrb.data(:,icol.sgrb.Epk);
    synSam.sgrb.T90 = sgrb.data(:,icol.sgrb.T90);
    synSam.sgrb.detProb = sgrb.data(:,icol.sgrb.detProb);
    
    % Create T90 bins
    synSam.all.binEdges = getBinEdges([synSam.sgrb.T90; synSam.lgrb.T90], nbins);
    synSam.all.binCenters = getBinCenters(synSam.all.binEdges);
    
    % Determine BATSE-detectable events and reduce dataset size to be comparable to BATSE LGRB catalog. Repeat 'nruns' times.
    for i = 1:nruns
        % If detection probability is above a uniformly generated random number, detect event
        Mask = synSam.lgrb.detProb > unifrnd(0, 1, length(synSam.lgrb.detProb), 1);
        synSam.lgrb.all.detectedVec(i).Epk = synSam.lgrb.Epk(Mask, 1);
        synSam.lgrb.all.detectedVec(i).T90 = synSam.lgrb.T90(Mask, 1);
        synSam.lgrb.all.detectedVec(i).Epk = synSam.lgrb.all.detectedVec(i).Epk(1:skip.lgrb:end);
        synSam.lgrb.all.detectedVec(i).T90 = synSam.lgrb.all.detectedVec(i).T90(1:skip.lgrb:end);
        Mask = synSam.sgrb.detProb > unifrnd(0, 1, length(synSam.sgrb.detProb), 1);
        synSam.sgrb.all.detectedVec(i).Epk = synSam.sgrb.Epk(Mask, 1);
        synSam.sgrb.all.detectedVec(i).T90 = synSam.sgrb.T90(Mask, 1);
        synSam.sgrb.all.detectedVec(i).Epk = synSam.sgrb.all.detectedVec(i).Epk(1:skip.sgrb:end);
        synSam.sgrb.all.detectedVec(i).T90 = synSam.sgrb.all.detectedVec(i).T90(1:skip.sgrb:end);
        synSam.both.all.detectedVec(i).Epk = [synSam.sgrb.all.detectedVec(i).Epk.', synSam.lgrb.all.detectedVec(i).Epk.'].';
        synSam.both.all.detectedVec(i).T90 = [synSam.sgrb.all.detectedVec(i).T90.', synSam.lgrb.all.detectedVec(i).T90.'].';
        synSam.both.all.detectedVec(i).countsVec = histcounts(synSam.both.all.detectedVec(i).T90, synSam.all.binEdges);
        
        % Cut the data set at the mean
        synSam.both.cut50.detectedVec(i).thresh = exp(mean(log(synSam.both.all.detectedVec(i).Epk(synSam.both.all.detectedVec(i).T90 >= 20))));
        Mask = synSam.both.all.detectedVec(i).Epk < synSam.both.cut50.detectedVec(i).thresh;
        synSam.both.cut50.detectedVec(i).Epk = synSam.both.all.detectedVec(i).Epk(Mask, 1);
        synSam.both.cut50.detectedVec(i).T90 = synSam.both.all.detectedVec(i).T90(Mask, 1);
        synSam.both.cut50.detectedVec(i).countsVec = histcounts(synSam.both.cut50.detectedVec(i).T90, synSam.all.binEdges);
        
        % Cut the data set one std above mean
        synSam.both.cut84.detectedVec(i).thresh = exp(mean(log(synSam.both.all.detectedVec(i).Epk(synSam.both.all.detectedVec(i).T90 >= 20))) ...
                                                    + std(log(synSam.both.all.detectedVec(i).Epk(synSam.both.all.detectedVec(i).T90 >= 20))));
        Mask = synSam.both.all.detectedVec(i).Epk < synSam.both.cut84.detectedVec(i).thresh;
        synSam.both.cut84.detectedVec(i).Epk = synSam.both.all.detectedVec(i).Epk(Mask, 1);
        synSam.both.cut84.detectedVec(i).T90 = synSam.both.all.detectedVec(i).T90(Mask, 1);
        synSam.both.cut84.detectedVec(i).countsVec = histcounts(synSam.both.cut84.detectedVec(i).T90, synSam.all.binEdges);
        
        % Cut the data set one std below mean
        synSam.both.cut16.detectedVec(i).thresh = exp(mean(log(synSam.both.all.detectedVec(i).Epk(synSam.both.all.detectedVec(i).T90 >= 20))) ...
                                                    - std(log(synSam.both.all.detectedVec(i).Epk(synSam.both.all.detectedVec(i).T90 >= 20))));
        Mask = synSam.both.all.detectedVec(i).Epk < synSam.both.cut16.detectedVec(i).thresh;
        synSam.both.cut16.detectedVec(i).Epk = synSam.both.all.detectedVec(i).Epk(Mask, 1);
        synSam.both.cut16.detectedVec(i).T90 = synSam.both.all.detectedVec(i).T90(Mask, 1);
        synSam.both.cut16.detectedVec(i).countsVec = histcounts(synSam.both.cut16.detectedVec(i).T90, synSam.all.binEdges);
    end
    clear('lgrb', 'sgrb', 'icol', 'Mask');
    synSam = rmfield(synSam, {'lgrb', 'sgrb'});
    
    % Grab a random set of detected data for Figure 1
    i = randi(nruns);
    synSam.both.all.detectedRand.Epk = synSam.both.all.detectedVec(i).Epk;
    synSam.both.all.detectedRand.T90 = synSam.both.all.detectedVec(i).T90;
    synSam.both.cut50.detectedRand.Epk = synSam.both.cut50.detectedVec(i).Epk;
    synSam.both.cut50.detectedRand.T90 = synSam.both.cut50.detectedVec(i).T90;
    synSam.both.cut84.detectedRand.Epk = synSam.both.cut84.detectedVec(i).Epk;
    synSam.both.cut84.detectedRand.T90 = synSam.both.cut84.detectedVec(i).T90;
    synSam.both.cut16.detectedRand.Epk = synSam.both.cut16.detectedVec(i).Epk;
    synSam.both.cut16.detectedRand.T90 = synSam.both.cut16.detectedVec(i).T90;
    synSam.both.cut50.detectedRand.thresh = synSam.both.cut50.detectedVec(i).thresh;
    synSam.both.cut84.detectedRand.thresh = synSam.both.cut84.detectedVec(i).thresh;
    synSam.both.cut16.detectedRand.thresh = synSam.both.cut16.detectedVec(i).thresh;
    
    for i = 1:nbins
        % Switch the positions of countsVec and detectedVec in the struct
        for j = 1:nruns
            synSam.both.all.countsVec(i).detectedVec(j) = synSam.both.all.detectedVec(j).countsVec(i);
            synSam.both.cut50.countsVec(i).detectedVec(j) = synSam.both.cut50.detectedVec(j).countsVec(i);
            synSam.both.cut84.countsVec(i).detectedVec(j) = synSam.both.cut84.detectedVec(j).countsVec(i);
            synSam.both.cut16.countsVec(i).detectedVec(j) = synSam.both.cut16.detectedVec(j).countsVec(i);
        end
        % Calculate the mean, the 95th percentile, and the 5th percentile
        synSam.both.all.percentile50.countsVec(i) = round(prctile(synSam.both.all.countsVec(i).detectedVec, 50));
        synSam.both.cut50.percentile50.countsVec(i) = round(prctile(synSam.both.cut50.countsVec(i).detectedVec, 50));
        synSam.both.cut84.percentile50.countsVec(i) = round(prctile(synSam.both.cut84.countsVec(i).detectedVec, 50));
        synSam.both.cut16.percentile50.countsVec(i) = round(prctile(synSam.both.cut16.countsVec(i).detectedVec, 50));
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
    synSam.both.all = rmfield(synSam.both.all, {'countsVec', 'detectedVec'});
    synSam.both.cut50 = rmfield(synSam.both.cut50, {'countsVec', 'detectedVec'});
    synSam.both.cut84 = rmfield(synSam.both.cut84, {'countsVec', 'detectedVec'});
    synSam.both.cut16 = rmfield(synSam.both.cut16, {'countsVec', 'detectedVec'});
    
    % Calculate ndata
    synSam.both.all.percentile50.ndata = sum(synSam.both.all.percentile50.countsVec);
    synSam.both.cut50.percentile50.ndata = sum(synSam.both.cut50.percentile50.countsVec);
    synSam.both.cut84.percentile50.ndata = sum(synSam.both.cut84.percentile50.countsVec);
    synSam.both.cut16.percentile50.ndata = sum(synSam.both.cut16.percentile50.countsVec);
    
    % Trim Zeros off of binned data for proper algorithmic x-axis limits
    [synSam.both.all.trim.percentile50.countsVec, synSam.both.all.trim.percentile50.binCenters, synSam.both.all.trim.percentile50.binEdges] ...
        = trimZeros(synSam.both.all.percentile50.countsVec, synSam.all.binEdges);
    
    % Convert the dependent variable dN/dlog(T) to dN/dT by dividing it by T90
    synSam.both.all.trim.percentile50.dNdT = synSam.both.all.trim.percentile50.countsVec ./ synSam.both.all.trim.percentile50.binCenters;
    synSam.both.cut50.percentile50.dNdT = synSam.both.cut50.percentile50.countsVec ./ synSam.all.binCenters;
    synSam.both.cut84.percentile50.dNdT = synSam.both.cut84.percentile50.countsVec ./ synSam.all.binCenters;
    synSam.both.cut16.percentile50.dNdT = synSam.both.cut16.percentile50.countsVec ./ synSam.all.binCenters;
    
    % Save mat file
    synSam.output.path = "../out"; if ~isfolder(synSam.output.path); mkdir(synSam.output.path); end
    save(synSam.output.path + "/" + matFileName, "synSam");
else
    % Load mat file
    synSam.output.path = "../out";
    load(synSam.output.path + "/" + matFileName); % loads synSam object   
end

% Figure 1: Epk vs T90 of random detected dataset
figure("color", figureColor); hold on; box on;
    scatter(synSam.both.all.detectedRand.T90, synSam.both.all.detectedRand.Epk, 5, defaultColor, 'filled', 'MarkerFaceAlpha', 0.25);
    scatter(synSam.both.cut84.detectedRand.T90, synSam.both.cut84.detectedRand.Epk, 5, defaultColor, 'filled', 'MarkerFaceAlpha', 0.25);
    scatter(synSam.both.cut50.detectedRand.T90, synSam.both.cut50.detectedRand.Epk, 5, defaultColor, 'filled', 'MarkerFaceAlpha', 0.25);
    scatter(synSam.both.cut16.detectedRand.T90, synSam.both.cut16.detectedRand.Epk, 5, defaultColor, 'filled');
    yline([synSam.both.cut84.detectedRand.thresh, synSam.both.cut50.detectedRand.thresh, synSam.both.cut16.detectedRand.thresh]);
    [xlower, xupper] = plotLimits(synSam.both.all.detectedRand.T90, 'log');
    xlim([xlower, xupper]);
    [ylower, yupper] = plotLimits(synSam.both.all.detectedRand.Epk, 'log');
    ylim([ylower, yupper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("E_{peak} [keV]", "interpreter", "tex", "fontsize", fontSize);
    title("Detectable GRBs (n = " + synSam.both.all.percentile50.ndata + ")", "fontSize", fontSize+1);
hold off;

% Figure 2: Plot dN/dT of full detected dataset
figure("color", figureColor); hold on; box on;
    hist2stairs(synSam.both.all.trim.percentile50.dNdT, synSam.both.all.trim.percentile50.binEdges, defaultColor, lineWidth);
    if mergeBins
        [newBins, binEdges, dNdT] = binMergeStairPlot( synSam.both.all.trim.percentile50.countsVec, synSam.both.all.trim.percentile50.binEdges ...
                                                     , secondaryColor, lineWidth);
        legend([synSam.both.all.percentile50.ndata + " detected GRBs" + newline + "spanning " + nbins + " bins" ...
               , "Adjacent bins with <5 events merged," + newline + "resulting in " + newBins + " bins"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Full Detected Data (mean of " + nruns + " runs)", "fontSize", fontSize+1);
        [xlower, xupper] = plotLimits(binEdges, 'log');
        xlim([xlower, xupper]);
        [ylower, yupper] = plotLimits(dNdT, 'log');
        ylim([ylower, yupper]);
    else
        legend([synSam.both.all.percentile50.ndata + " detected GRBs," + newline + "mean of " + nruns + " runs"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Full Detected Data", "fontSize", fontSize+1);
        [xlower, xupper] = plotLimits(synSam.both.all.trim.percentile50.binEdges, 'log');
        xlim([xlower, xupper]);
        [ylower, yupper] = plotLimits(synSam.both.all.trim.percentile50.dNdT, 'log');
        ylim([ylower, yupper]);
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
hold off;

% Figure 3: Plot dN/dT of dataset cut at mean
figure("color", figureColor); hold on; box on;
    hist2stairs(synSam.both.cut50.percentile50.dNdT, synSam.all.binEdges, defaultColor, lineWidth);
    if mergeBins
        [newBins, ~, dNdT] = binMergeStairPlot(synSam.both.cut50.percentile50.countsVec, synSam.all.binEdges, secondaryColor, lineWidth);
        legend([synSam.both.cut50.percentile50.ndata + " detected GRBs" + newline + "spanning " + nbins + " bins" ...
               , "Adjacent bins with <5 events merged," + newline + "resulting in " + newBins + " bins"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean (mean of " + nruns + " runs)", "fontSize", fontSize+1);
        [ylower, yupper] = plotLimits(dNdT, 'log');
        ylim([ylower, yupper]);
    else
        legend([synSam.both.cut50.percentile50.ndata + " detected GRBs," + newline + "mean of " + nruns + " runs"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean (for T90 >= 20 s)", "fontSize", fontSize+1);
        [ylower, yupper] = plotLimits(synSam.both.cut50.percentile50.dNdT, 'log');
        ylim([ylower, yupper]);
    end
    xlim([xlower, xupper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
hold off;

% Figure 4: Plot dN/dT of dataset cut at mean + std
figure("color", figureColor); hold on; box on;
    hist2stairs(synSam.both.cut84.percentile50.dNdT, synSam.all.binEdges, defaultColor, lineWidth);
    if mergeBins
        [newBins, ~, dNdT] = binMergeStairPlot(synSam.both.cut84.percentile50.countsVec, synSam.all.binEdges, secondaryColor, lineWidth);
        legend([synSam.both.cut84.percentile50.ndata + " detected GRBs" + newline + "spanning " + nbins + " bins" ...
               , "Adjacent bins with <5 events merged," + newline + "resulting in " + newBins + " bins"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean + Std (mean of " + nruns + " runs)", "fontSize", fontSize+1);
        [ylower, yupper] = plotLimits(dNdT, 'log');
        ylim([ylower, yupper]);
    else
        legend([synSam.both.cut84.percentile50.ndata + " detected GRBs," + newline + "mean of " + nruns + " runs"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean + Std (for T90 >= 20 s)", "fontSize", fontSize+1);
        [ylower, yupper] = plotLimits(synSam.both.cut84.percentile50.dNdT, 'log');
        ylim([ylower, yupper]);
    end
    xlim([xlower, xupper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
hold off;

% Figure 5: Plot dN/dT of dataset cut at mean - std
figure("color", figureColor); hold on; box on;
    hist2stairs(synSam.both.cut16.percentile50.dNdT, synSam.all.binEdges, defaultColor, lineWidth);
    if mergeBins
        [newBins, ~, dNdT] = binMergeStairPlot(synSam.both.cut16.percentile50.countsVec, synSam.all.binEdges, secondaryColor, lineWidth);
        legend([synSam.both.cut16.percentile50.ndata + " detected GRBs" + newline + "spanning " + nbins + " bins" ...
               , "Adjacent bins with <5 events merged," + newline + "resulting in " + newBins + " bins"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean - Std (mean of " + nruns + " runs)", "fontSize", fontSize+1);
        [ylower, yupper] = plotLimits(dNdT, 'log');
        ylim([ylower, yupper]);
    else
        legend([synSam.both.cut16.percentile50.ndata + " detected GRBs," + newline + "mean of " + nruns + " runs"] ...
               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        title("Detected Data Cut at Epk Mean - Std (for T90 >= 20 s)", "fontSize", fontSize+1);
        [ylower, yupper] = plotLimits(synSam.both.cut16.percentile50.dNdT, 'log');
        ylim([ylower, yupper]);
    end
    xlim([xlower, xupper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
hold off;