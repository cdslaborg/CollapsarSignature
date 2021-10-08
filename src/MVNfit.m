% clear all;
close all;
% filePath = mfilename('fullpath');
% [currentDir,fileName,fileExt] = fileparts(filePath); cd(currentDir);
% cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
% addpath([genpath("../../../libmatlab"), genpath("./functions")], "-begin");

% Figure Parameters
fontSize = 14;
lineWidth = 1.2;
figureColor = "white";
defaultColor = [0, 0.4470, 0.7410];
secondaryColor = 'green'; %'#F80'; %orange
sgrbColor = 'blue';
lgrbColor = 'red';

% Options
nsample = 2000; % number of samples drawn from PDF
nbins = 45; % number of bins in which to bin T90 data
mergeBins = true; % merge adjacent bins with <5 events
saveNewImages = false; % export figure on or off. Must have export_fig installed.

% BATSE T90 & Epk Data
% batse = importdata("..\in\Batse_orig.xlsx");
% icol.T90 = 14; % column index of T90
% icol.log10Epk = 17;
% data.batse.T90 = batse.data.x1966GRBs(:, icol.T90);
% data.batse.Epk = 10.^(batse.data.x1966GRBs(:, icol.log10Epk));
% data.batse.logVec = [log(data.batse.T90), log(data.batse.Epk)];
% clear('batse', 'icol');

% Fit Double MVN - Must have Statistics and Machine Learning Toolbox installed
% GMModel = fitgmdist(data.batse.logVec, 2); % specify two-component MVN
gmPDF = @(x, y) arrayfun(@(x0, y0) pdf(GMModel, [x0 y0]), log(x), log(y));

% Sample the PDF
% data.synSam.all.sgrbVec = exp(mvnrnd(GMModel.mu(1,:), GMModel.Sigma(:,:,1), round(GMModel.ComponentProportion(1)*nsample)));
% data.synSam.all.lgrbVec = exp(mvnrnd(GMModel.mu(2,:), GMModel.Sigma(:,:,2), round(GMModel.ComponentProportion(2)*nsample)));
% data.synSam.all.bothVec = [data.synSam.all.sgrbVec(:,:); data.synSam.all.lgrbVec(:,:)];

% Cut Data
% criteria = data.synSam.all.bothVec(:,1) >= 20; % only consider where T90 >= 20 s
% data.synSam.cut.thresh = exp(mean(log(data.synSam.all.bothVec(criteria, 2))));
% criteria = data.synSam.all.sgrbVec(:,2) < data.synSam.cut.thresh; % only consider where Epk < thresh
% data.synSam.cut.sgrbVec = data.synSam.all.sgrbVec(criteria, :);
% criteria = data.synSam.all.lgrbVec(:,2) < data.synSam.cut.thresh; % only consider where Epk < thresh
% data.synSam.cut.lgrbVec = data.synSam.all.lgrbVec(criteria, :);
% data.synSam.cut.bothVec = [data.synSam.cut.sgrbVec(:,:); data.synSam.cut.lgrbVec(:,:)];
% clear('criteria');

% Bin Data
data.synSam.all.T90.binEdges = getBinEdges(data.synSam.all.bothVec(:,1), nbins);
data.synSam.all.T90.binWidths = getBinWidths(data.synSam.all.T90.binEdges);
data.synSam.all.T90.countsVec = histcounts(data.synSam.all.bothVec(:,1), data.synSam.all.T90.binEdges);
data.synSam.cut.T90.countsVec = histcounts(data.synSam.cut.bothVec(:,1), data.synSam.all.T90.binEdges);

% Trim Zeros off of binned data for proper algorithmic x-axis limits
[data.synSam.all.T90.trim.countsVec, data.synSam.all.T90.trim.binWidths, data.synSam.all.T90.trim.binEdges] ...
    = trimZeros(data.synSam.all.T90.countsVec, data.synSam.all.T90.binEdges);
[data.synSam.cut.T90.trim.countsVec, data.synSam.cut.T90.trim.binWidths, data.synSam.cut.T90.trim.binEdges] ...
    = trimZeros(data.synSam.cut.T90.countsVec, data.synSam.all.T90.binEdges);

% Convert the histogram counts to dN/dT by dividing each bin by its width
data.synSam.all.T90.trim.dNdT = data.synSam.all.T90.trim.countsVec ./ data.synSam.all.T90.trim.binWidths;
data.synSam.cut.T90.trim.dNdT = data.synSam.cut.T90.trim.countsVec ./ data.synSam.cut.T90.trim.binWidths;

% Figure 1: BATSE Epk vs. T90
figure("color", figureColor); hold on; box on;
    scatter(data.batse.T90, data.batse.Epk, 5, defaultColor, 'filled')
    g = gca;
        set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
        xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
        ylabel("E_{peak} [keV]", "interpreter", "tex", "fontsize", fontSize);
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(data.batse.T90, 'log');
        xlim([ax.lower, ax.upper]);
        xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(data.batse.Epk, 'log');
        ylim([ax.lower, ax.upper]);
        yticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    fcontour(gmPDF, [g.XLim g.YLim])
hold off;

% Figure 2: Synthetic Data
figure("color", figureColor); hold on; box on;
    scatter(data.synSam.cut.sgrbVec(:,1), data.synSam.cut.sgrbVec(:,2), 4, sgrbColor, 'filled');
    scatter(data.synSam.cut.lgrbVec(:,1), data.synSam.cut.lgrbVec(:,2), 4, lgrbColor, 'filled');
    yline(data.synSam.cut.thresh, "color", "black", "LineWidth", lineWidth);
    scatter(data.synSam.all.sgrbVec(:,1), data.synSam.all.sgrbVec(:,2), 4, sgrbColor, 'filled', 'MarkerFaceAlpha', 0.4);
    scatter(data.synSam.all.lgrbVec(:,1), data.synSam.all.lgrbVec(:,2), 4, lgrbColor, 'filled', 'MarkerFaceAlpha', 0.4);
    leg = legend(["SGRBs (n = " + length(data.synSam.all.sgrbVec) + ")", "LGRBs (n = " + length(data.synSam.all.lgrbVec) + ")" ...
                 , "E_{pk} threshold = " + round(data.synSam.cut.thresh) + " keV"], "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
        leg.BoxFace.ColorType='truecoloralpha';
        leg.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
        xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
        ylabel("E_{peak} [keV]", "interpreter", "tex", "fontsize", fontSize);
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(data.batse.T90, 'log'); %[data.synSam.all.sgrbVec(:,1); data.synSam.all.lgrbVec(:,1)]
        xlim([ax.lower, ax.upper]);
        xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(data.batse.Epk, 'log'); %[data.synSam.all.sgrbVec(:,2); data.synSam.all.lgrbVec(:,2)]
        ylim([ax.lower, ax.upper]);
        yticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
hold off;

% Figure 3: Synthetic T90 Histograms
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', data.synSam.all.T90.binEdges, 'BinCounts', data.synSam.all.T90.countsVec, 'FaceColor', defaultColor);
    histogram('BinEdges', data.synSam.all.T90.binEdges, 'BinCounts', data.synSam.cut.T90.countsVec, 'FaceColor', secondaryColor);
    legend(["Synthetic T_{90} (all)", "Synthetic T_{90} (E_{pk} < threshold)"], "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
        xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
        ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(data.synSam.all.T90.binEdges, 'log');
        xlim([ax.lower, ax.upper]);
        xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
        [~, ax.upper, ~, ~, ~] = plotLimits(data.synSam.all.T90.countsVec, 'linear');
        ylim([0, ax.upper]);
hold off;

% Figure 4: Synthetic dN/DT
figure("color", figureColor); hold on; box on;
    if mergeBins
        [nbinsNew1, binEdges1, dNdT1] = binMergeStairPlot(data.synSam.all.T90.trim.countsVec, data.synSam.all.T90.trim.binEdges, defaultColor, lineWidth);
        [nbinsNew2, binEdges2, dNdT2] = binMergeStairPlot(data.synSam.cut.T90.trim.countsVec, data.synSam.cut.T90.trim.binEdges, secondaryColor, lineWidth);
        [ax.xlower, ax.xupper, ax.xlowTick, ax.xupTick, ax.xspan] = plotLimits([binEdges1, binEdges2], 'log');
        [ax.ylower, ax.yupper, ax.ylowTick, ax.yupTick, ax.yspan] = plotLimits([dNdT1, dNdT2], 'log');
        textpos.x = exp(log(ax.xupper) - (log(ax.xupper) - log(ax.xlower))*0.03);
        textpos.y = exp(log(ax.yupper) - (log(ax.yupper) - log(ax.ylower))*0.03);
        text(textpos.x, textpos.y, 'Bin Merging On', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontSize', fontSize-3);
    else
        hist2stairs(data.synSam.all.T90.trim.dNdT, data.synSam.all.T90.trim.binEdges, defaultColor, lineWidth);
        hist2stairs(data.synSam.cut.T90.trim.dNdT, data.synSam.cut.T90.trim.binEdges, secondaryColor, lineWidth);
        [ax.xlower, ax.xupper, ax.xlowTick, ax.xupTick, ax.xspan] = plotLimits([data.synSam.all.T90.trim.binEdges, data.synSam.cut.T90.trim.binEdges], 'log');
        [ax.ylower, ax.yupper, ax.ylowTick, ax.yupTick, ax.yspan] = plotLimits([data.synSam.all.T90.trim.dNdT, data.synSam.cut.T90.trim.dNdT], 'log');
        textpos.x = exp(log(ax.xupper) - (log(ax.xupper) - log(ax.xlower))*0.03);
        textpos.y = exp(log(ax.yupper) - (log(ax.yupper) - log(ax.ylower))*0.03);
        text(textpos.x, textpos.y, 'Bin Merging Off', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontSize', fontSize-3);
    end
    legend(["Synthetic T_{90} (all)", "Synthetic T_{90} (E_{pk} < threshold)"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
        xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
        ylabel("dN / dT_{90}", "interpreter", "tex", "fontsize", fontSize);
        xlim([ax.xlower, ax.xupper]);
        xticks(10.^(linspace(ax.xlowTick, ax.xupTick, ax.xspan)));
        ylim([ax.ylower, ax.yupper]);
        yticks(10.^(linspace(ax.ylowTick, ax.yupTick, ax.yspan)));
hold off;