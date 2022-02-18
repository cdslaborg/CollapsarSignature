clear all;
close all;
filePath = mfilename('fullpath');
[currentDir,fileName,fileExt] = fileparts(filePath); cd(currentDir);
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath([genpath("../../../libmatlab"), genpath("./functions")], "-begin");

% Figure Parameters
fontSize = 14;
lineWidth = 1.2;
figureColor = "white";
defaultColor = [0, 0.4470, 0.7410];
secondaryColor = '#F80'; %orange
sgrbColor = 'blue';
lgrbColor = 'red';
legendAlpha = 0.7;

% Options
xvar = 'T90'; % choose x-axis variable. Must be 'T90' or 'T50'.
yvar = 'HR32'; % choose y-axis variable. Must be 'HR32', 'HR4321', or 'Epk'.
threshold = 5; % value of (T90 / T50) to cut the dataset, discarding values below this threshold as not likely LGRBs. Default is 5 seconds.
facInc = 0; % factor by which to increase datset size for better fitting. Default is 0. Must be an integer.
nsample = 2000; % number of samples drawn from PDF
nbins = 50; % number of bins in which to bin xvar data
mergeBins = true; % merge adjacent bins with <5 events
saveNewImages = true; % export figure on or off. Must have export_fig installed.

if ~any(strcmpi(xvar, ["T90", "T50"])); error("'xvar' is not an allowable value"); end
if ~any(strcmpi(yvar, ["HR32", "HR4321", "Epk"])); error("'yvar' is not an allowable value"); end

% BATSE T90 & Epk Data
batse = importdata("..\in\Batse_orig.xlsx");
icol.f1 = 2; % column index of Fluence Channel 1
icol.f2 = 4;
icol.f3 = 6;
icol.f4 = 8;
icol.T50 = 11;
icol.T90 = 14;
icol.log10Epk = 17;
f1.values = batse.data.x1966GRBs(:, icol.f1);
f2.values = batse.data.x1966GRBs(:, icol.f2);
f3.values = batse.data.x1966GRBs(:, icol.f3);
f4.values = batse.data.x1966GRBs(:, icol.f4);
T50.values = batse.data.x1966GRBs(:, icol.T50);
T90.values = batse.data.x1966GRBs(:, icol.T90);
Epk.values = 10.^(batse.data.x1966GRBs(:, icol.log10Epk));

% correctly set 'data.batse.xvar' and 'data.batse.yvar'
data.batse.labels = {"x", "y"}; %#ok<*CLARRSTR> 
if any(strcmpi(yvar, ["HR32", "HR4321"]))
    f1.nonzeroInd = f1.values ~= 0.;
    f2.nonzeroInd = f2.values ~= 0.;
    f3.nonzeroInd = f3.values ~= 0.;
    f4.nonzeroInd = f4.values ~= 0.;
    if strcmpi(yvar, "HR32")
        nonzeroInd = f2.nonzeroInd & f3.nonzeroInd;
        data.batse.yvar = f3.values(nonzeroInd) ./ f2.values(nonzeroInd);
        data.batse.labels{2} = "HR_{32}";
    else
        nonzeroInd = f1.nonzeroInd & f2.nonzeroInd & f3.nonzeroInd & f4.nonzeroInd;
        data.batse.yvar = (f4.values(nonzeroInd) + f3.values(nonzeroInd)) ./ (f2.values(nonzeroInd) + f1.values(nonzeroInd));
        data.batse.labels{2} = "HR_{4321}";
    end
else
    data.batse.yvar = Epk.values;
    data.batse.labels{2} = "E_{pk} [keV]";
end
if strcmpi(xvar, "T90")
    data.batse.labels{1} = "T_{90} [s]";
    if exist("nonzeroInd", "var")
        data.batse.xvar = T90.values(nonzeroInd);
    else
        data.batse.xvar = T90.values;
    end
else
    data.batse.labels{1} = "T_{50} [s]";
    if exist("nonzeroInd", "var")
        data.batse.xvar = T50.values(nonzeroInd);
    else
        data.batse.xvar = T50.values;
    end
end
data.batse.all.ndata = length(data.batse.xvar);
clear('batse', 'icol', 'f1', 'f2', 'f3', 'f4', 'T50', 'T90', 'Epk', 'nonzeroInd');

% arrange data in a vector for MVN fitting
data.batse.all.vec = [data.batse.xvar, data.batse.yvar];
data.batse.all.logVec = log(data.batse.all.vec);
data.batse = rmfield(data.batse, {'xvar', 'yvar'});

% Increase the dataset size for better MVN fitting
if facInc ~= 0
    percent = 0.01;
    for i = 1:facInc
        varVec = normrnd(0, 0.5, [data.batse.all.ndata, 2]);
        for col = 1:2
            logMin = min(data.batse.all.logVec(:,col));
            logMax = max(data.batse.all.logVec(:,col));
            logRange = logMax - logMin;
            vary = logRange * percent;
            varVec(:,col) = varVec(:,col) * vary;
        end
        varVec = data.batse.all.logVec(data.batse.all.ndata*(i-1)+1:data.batse.all.ndata*i,:) + varVec;
        data.batse.all.logVec = [data.batse.all.logVec; varVec];
    end
    data.batse.all.vec = exp(data.batse.all.logVec);
    clear('percent', 'varVec', 'col', 'logMin', 'logMax', 'logRange', 'vary');
end

% Fit Double MVN - Must have Statistics and Machine Learning Toolbox installed
if strcmpi(xvar, 'T90') && strcmpi(yvar, 'HR32')
    load("../out/modelParamHR32T90.mat");
    S.mu = GMModel.mu;
    S.Sigma = GMModel.Sigma;
    S.ComponentProportion = GMModel.ComponentProportion;
    clear('GMModel');
    options = statset('MaxIter',10000); ... 'Display','final',
    GMModel = fitgmdist(data.batse.all.logVec, 2, 'Options', options, 'Replicates', 100); ...'Start', S);
else
    GMModel = fitgmdist(data.batse.all.logVec, 2); % specify two-component MVN
end
gmPDF = @(x, y) arrayfun(@(x0, y0) pdf(GMModel, [x0 y0]), log(x), log(y));

% Sample the PDF
dummy1 = exp(mvnrnd(GMModel.mu(1,:), GMModel.Sigma(:,:,1), round(GMModel.ComponentProportion(1)*nsample)));
dummy2 = exp(mvnrnd(GMModel.mu(2,:), GMModel.Sigma(:,:,2), round(GMModel.ComponentProportion(2)*nsample)));
if mean(dummy1(:, 1)) < mean(dummy2(:, 1))
    data.synSam.all.sgrbVec = dummy1;
    data.synSam.all.lgrbVec = dummy2;
else
    data.synSam.all.sgrbVec = dummy2;
    data.synSam.all.lgrbVec = dummy1;
end
data.synSam.all.bothVec = [data.synSam.all.sgrbVec; data.synSam.all.lgrbVec];
clear('dummy1', 'dummy2');

% Cut Data: only consider where xvar >= threshold
criteria = data.batse.all.vec(:, 1) >= threshold;
data.batse.cutY.thresh = exp(median(log(data.batse.all.vec(criteria, 2))));
criteria = data.synSam.all.bothVec(:, 1) >= threshold;
data.synSam.cutY.thresh = exp(median(log(data.synSam.all.bothVec(criteria, 2))));

% Cut Data: only consider where yvar < thresh
criteria = data.batse.all.vec(:, 2) < data.batse.cutY.thresh;
data.batse.cutY.vec = data.batse.all.vec(criteria, :);
criteria = data.synSam.all.sgrbVec(:, 2) < data.synSam.cutY.thresh;
data.synSam.cutY.sgrbVec = data.synSam.all.sgrbVec(criteria, :);
criteria = data.synSam.all.lgrbVec(:, 2) < data.synSam.cutY.thresh;
data.synSam.cutY.lgrbVec = data.synSam.all.lgrbVec(criteria, :);
data.synSam.cutY.bothVec = [data.synSam.cutY.sgrbVec; data.synSam.cutY.lgrbVec];
clear('criteria');

% Bin Data
data.batse.all.xvar.binEdges = getBinEdges(data.batse.all.vec(:, 1), nbins);
data.batse.all.xvar.binWidths = getBinWidths(data.batse.all.xvar.binEdges);
data.batse.all.xvar.countsVec = histcounts(data.batse.all.vec(:, 1), data.batse.all.xvar.binEdges);
data.batse.cutY.xvar.countsVec = histcounts(data.batse.cutY.vec(:, 1), data.batse.all.xvar.binEdges);
data.synSam.all.xvar.binEdges = getBinEdges(data.synSam.all.bothVec(:, 1), nbins);
data.synSam.all.xvar.binWidths = getBinWidths(data.synSam.all.xvar.binEdges);
data.synSam.all.xvar.countsVec = histcounts(data.synSam.all.bothVec(:, 1), data.synSam.all.xvar.binEdges);
data.synSam.cutY.xvar.countsVec = histcounts(data.synSam.cutY.bothVec(:, 1), data.synSam.all.xvar.binEdges);

% Trim Zeros off of binned data for proper algorithmic x-axis limits
[data.batse.all.xvar.trim.countsVec, data.batse.all.xvar.trim.binWidths, data.batse.all.xvar.trim.binEdges] ...
    = trimZeros(data.batse.all.xvar.countsVec, data.batse.all.xvar.binEdges);
[data.batse.cutY.xvar.trim.countsVec, data.batse.cutY.xvar.trim.binWidths, data.batse.cutY.xvar.trim.binEdges] ...
    = trimZeros(data.batse.cutY.xvar.countsVec, data.batse.all.xvar.binEdges);
[data.synSam.all.xvar.trim.countsVec, data.synSam.all.xvar.trim.binWidths, data.synSam.all.xvar.trim.binEdges] ...
    = trimZeros(data.synSam.all.xvar.countsVec, data.synSam.all.xvar.binEdges);
[data.synSam.cutY.xvar.trim.countsVec, data.synSam.cutY.xvar.trim.binWidths, data.synSam.cutY.xvar.trim.binEdges] ...
    = trimZeros(data.synSam.cutY.xvar.countsVec, data.synSam.all.xvar.binEdges);

% Convert the histogram counts to dN/dT by dividing each bin by its width
data.batse.all.xvar.trim.dNdT = data.batse.all.xvar.trim.countsVec ./ data.batse.all.xvar.trim.binWidths;
data.batse.cutY.xvar.trim.dNdT = data.batse.cutY.xvar.trim.countsVec ./ data.batse.cutY.xvar.trim.binWidths;
data.synSam.all.xvar.trim.dNdT = data.synSam.all.xvar.trim.countsVec ./ data.synSam.all.xvar.trim.binWidths;
data.synSam.cutY.xvar.trim.dNdT = data.synSam.cutY.xvar.trim.countsVec ./ data.synSam.cutY.xvar.trim.binWidths;

if facInc > 0
    colorVec = {'red', 'blue', 'green'};
else
    colorVec = {defaultColor};
end

% Figure 1: BATSE (HR32 / HR4321 / Epk) vs. (T90 / T50)
fig1 = figure("color", figureColor); hold on; box on;
    for i = 1:min([facInc + 1, 4])
        scatter( data.batse.all.vec(data.batse.all.ndata*(i-1)+1:data.batse.all.ndata*i, 1) ...
               , data.batse.all.vec(data.batse.all.ndata*(i-1)+1:data.batse.all.ndata*i, 2) ...
               , 4, colorVec{i}, 'filled')
    end
    g = gca;
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
        xlabel(data.batse.labels(1), "interpreter", "tex", "fontsize", fontSize);
        ylabel(data.batse.labels(2), "interpreter", "tex", "fontsize", fontSize);
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(data.batse.all.vec(:, 1), 'log');
        xlim([ax.lower, ax.upper]);
        xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(data.batse.all.vec(:, 2), 'log');
        ylim([ax.lower, ax.upper]);
        yticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    fcontour(gmPDF, [g.XLim g.YLim])
    if facInc == 0
        leg = legend(["BATSE GRBs (n = " + data.batse.all.ndata + ")", "Contour of MVN fit"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
            leg.BoxFace.ColorType='truecoloralpha';
            leg.BoxFace.ColorData=uint8(255*[1 1 1 legendAlpha]');
    end
    if saveNewImages
        export_fig("../out/Fig1_BATSE_" + yvar + "vs" + xvar + ".png", "-m4 -transparent");
    end
hold off;

% Figure 2: BATSE (T90 / T50) Histograms from (HR32 / HR4321 / Epk) cut
fig2 = figure("color", figureColor); hold on; box on;
    histogram('BinEdges', data.batse.all.xvar.binEdges, 'BinCounts', data.batse.all.xvar.countsVec, 'FaceColor', defaultColor);
    histogram('BinEdges', data.batse.all.xvar.binEdges, 'BinCounts', data.batse.cutY.xvar.countsVec, 'FaceColor', secondaryColor);
    leg = legend( ["BATSE " + data.batse.labels(1) + " (all)", "BATSE " + data.batse.labels(1) + " (" + data.batse.labels(2) + " < threshold)"] ...
          , "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
        leg.BoxFace.ColorType='truecoloralpha';
        leg.BoxFace.ColorData=uint8(255*[1 1 1 legendAlpha]');
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
        xlabel(data.batse.labels(1), "interpreter", "tex", "fontsize", fontSize);
        ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(data.batse.all.xvar.binEdges, 'log');
        xlim([ax.lower, ax.upper]);
        xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
        [~, ax.upper, ~, ~, ~] = plotLimits(data.batse.all.xvar.countsVec, 'linear');
        ylim([0, ax.upper]);
    if saveNewImages
        export_fig("../out/Fig2_BATSE_" + xvar + "Histograms.png", "-m4 -transparent");
    end
hold off;

% Figure 3: BATSE dN/DT from (HR32 / HR4321 / Epk) cut
fig3 = figure("color", figureColor); hold on; box on;
    if mergeBins
        [nbinsNew1, binEdges1, dNdT1] = binMergeStairPlot(data.batse.all.xvar.trim.countsVec, data.batse.all.xvar.trim.binEdges, defaultColor, lineWidth);
        [nbinsNew2, binEdges2, dNdT2] = binMergeStairPlot(data.batse.cutY.xvar.trim.countsVec, data.batse.cutY.xvar.trim.binEdges, secondaryColor, lineWidth);
        [ax.xlower, ax.xupper, ax.xlowTick, ax.xupTick, ax.xspan] = plotLimits([binEdges1, binEdges2], 'log');
        [ax.ylower, ax.yupper, ax.ylowTick, ax.yupTick, ax.yspan] = plotLimits([dNdT1, dNdT2], 'log');
        textpos.x = exp(log(ax.xupper) - (log(ax.xupper) - log(ax.xlower))*0.03);
        textpos.y = exp(log(ax.yupper) - (log(ax.yupper) - log(ax.ylower))*0.03);
        text(textpos.x, textpos.y, 'Bin Merging On', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontSize', fontSize-3);
    else
        hist2stairs(data.batse.all.xvar.trim.dNdT, data.batse.all.xvar.trim.binEdges, defaultColor, lineWidth);
        hist2stairs(data.batse.cutY.xvar.trim.dNdT, data.batse.cutY.xvar.trim.binEdges, secondaryColor, lineWidth);
        [ax.xlower, ax.xupper, ax.xlowTick, ax.xupTick, ax.xspan] = plotLimits([data.batse.all.xvar.trim.binEdges, data.batse.cutY.xvar.trim.binEdges], 'log');
        [ax.ylower, ax.yupper, ax.ylowTick, ax.yupTick, ax.yspan] = plotLimits([data.batse.all.xvar.trim.dNdT, data.batse.cutY.xvar.trim.dNdT], 'log');
        textpos.x = exp(log(ax.xupper) - (log(ax.xupper) - log(ax.xlower))*0.03);
        textpos.y = exp(log(ax.yupper) - (log(ax.yupper) - log(ax.ylower))*0.03);
        text(textpos.x, textpos.y, 'Bin Merging Off', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontSize', fontSize-3);
    end
    legend( ["BATSE " + data.batse.labels(1) + " (all)", "BATSE " + data.batse.labels(1) + " (" + data.batse.labels(2) + " < threshold)"] ...
          , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
        xlabel(data.batse.labels(1), "interpreter", "tex", "fontsize", fontSize);
        ylabel("dN / dT", "interpreter", "tex", "fontsize", fontSize);
        xlim([ax.xlower, ax.xupper]);
        xticks(10.^(linspace(ax.xlowTick, ax.xupTick, ax.xspan)));
        ylim([ax.ylower, ax.yupper]);
        yticks(10.^(linspace(ax.ylowTick, ax.yupTick, ax.yspan)));
    if saveNewImages
        export_fig("../out/Fig3_BATSE_dNdTvs" + xvar + ".png", "-m4 -transparent");
    end
hold off;

% Figure 4: Synthetic (HR32 / HR4321 / Epk) vs. (T90 / T50)
fig4 = figure("color", figureColor); hold on; box on;
    scatter(data.synSam.cutY.sgrbVec(:, 1), data.synSam.cutY.sgrbVec(:, 2), 4, sgrbColor, 'filled');
    scatter(data.synSam.cutY.lgrbVec(:, 1), data.synSam.cutY.lgrbVec(:, 2), 4, lgrbColor, 'filled');
    yline(data.synSam.cutY.thresh, "color", "black", "LineWidth", lineWidth);
    scatter(data.synSam.all.sgrbVec(:, 1), data.synSam.all.sgrbVec(:, 2), 4, sgrbColor, 'filled', 'MarkerFaceAlpha', 0.4);
    scatter(data.synSam.all.lgrbVec(:, 1), data.synSam.all.lgrbVec(:, 2), 4, lgrbColor, 'filled', 'MarkerFaceAlpha', 0.4);
    leg = legend( ["Synthetic SGRBs (n = " + length(data.synSam.all.sgrbVec) + ")", "Synthetic LGRBs (n = " + length(data.synSam.all.lgrbVec) + ")" ...
                , data.batse.labels(2) + " threshold = " + round(data.synSam.cutY.thresh, 1)], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        leg.BoxFace.ColorType='truecoloralpha';
        leg.BoxFace.ColorData=uint8(255*[1 1 1 legendAlpha]');
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
        xlabel(data.batse.labels(1), "interpreter", "tex", "fontsize", fontSize);
        ylabel(data.batse.labels(2), "interpreter", "tex", "fontsize", fontSize);
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(data.batse.all.vec(:, 1), 'log');
        xlim([ax.lower, ax.upper]);
        xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(data.batse.all.vec(:, 2), 'log');
        ylim([ax.lower, ax.upper]);
        yticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    if saveNewImages
        export_fig("../out/Fig4_SynSam_" + yvar + "vs" + xvar + ".png", "-m4 -transparent");
    end
hold off;

% Figure 5: Synthetic (T90 / T50) Histograms from (HR32 / HR4321 / Epk) cut
fig5 = figure("color", figureColor); hold on; box on;
    histogram('BinEdges', data.synSam.all.xvar.binEdges, 'BinCounts', data.synSam.all.xvar.countsVec, 'FaceColor', defaultColor);
    histogram('BinEdges', data.synSam.all.xvar.binEdges, 'BinCounts', data.synSam.cutY.xvar.countsVec, 'FaceColor', secondaryColor);
    leg = legend( ["Synthetic " + data.batse.labels(1) + " (all)", "Synthetic " + data.batse.labels(1) + " (" + data.batse.labels(2) + " < threshold)"] ...
          , "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
        leg.BoxFace.ColorType='truecoloralpha';
        leg.BoxFace.ColorData=uint8(255*[1 1 1 legendAlpha]');
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
        xlabel(data.batse.labels(1), "interpreter", "tex", "fontsize", fontSize);
        ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(data.synSam.all.xvar.binEdges, 'log');
        xlim([ax.lower, ax.upper]);
        xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
        [~, ax.upper, ~, ~, ~] = plotLimits(data.synSam.all.xvar.countsVec, 'linear');
        ylim([0, ax.upper]);
    if saveNewImages
        export_fig("../out/Fig5_SynSam_" + xvar + "Histograms.png", "-m4 -transparent");
    end
hold off;

% Figure 6: Synthetic dN/DT from (HR32 / HR4321 / Epk) cut
fig6 = figure("color", figureColor); hold on; box on;
    if mergeBins
        [nbinsNew1, binEdges1, dNdT1] = binMergeStairPlot(data.synSam.all.xvar.trim.countsVec, data.synSam.all.xvar.trim.binEdges, defaultColor, lineWidth);
        [nbinsNew2, binEdges2, dNdT2] = binMergeStairPlot(data.synSam.cutY.xvar.trim.countsVec, data.synSam.cutY.xvar.trim.binEdges, secondaryColor, lineWidth);
        [ax.xlower, ax.xupper, ax.xlowTick, ax.xupTick, ax.xspan] = plotLimits([binEdges1, binEdges2], 'log');
        [ax.ylower, ax.yupper, ax.ylowTick, ax.yupTick, ax.yspan] = plotLimits([dNdT1, dNdT2], 'log');
        textpos.x = exp(log(ax.xupper) - (log(ax.xupper) - log(ax.xlower))*0.03);
        textpos.y = exp(log(ax.yupper) - (log(ax.yupper) - log(ax.ylower))*0.03);
        text(textpos.x, textpos.y, 'Bin Merging On', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontSize', fontSize-3);
    else
        hist2stairs(data.synSam.all.xvar.trim.dNdT, data.synSam.all.xvar.trim.binEdges, defaultColor, lineWidth);
        hist2stairs(data.synSam.cutY.xvar.trim.dNdT, data.synSam.cutY.xvar.trim.binEdges, secondaryColor, lineWidth);
        [ax.xlower, ax.xupper, ax.xlowTick, ax.xupTick, ax.xspan] = plotLimits([data.synSam.all.xvar.trim.binEdges, data.synSam.cutY.xvar.trim.binEdges], 'log');
        [ax.ylower, ax.yupper, ax.ylowTick, ax.yupTick, ax.yspan] = plotLimits([data.synSam.all.xvar.trim.dNdT, data.synSam.cutY.xvar.trim.dNdT], 'log');
        textpos.x = exp(log(ax.xupper) - (log(ax.xupper) - log(ax.xlower))*0.03);
        textpos.y = exp(log(ax.yupper) - (log(ax.yupper) - log(ax.ylower))*0.03);
        text(textpos.x, textpos.y, 'Bin Merging Off', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'fontSize', fontSize-3);
    end
    legend( ["Synthetic " + data.batse.labels(1) + " (all)", "Synthetic " + data.batse.labels(1) + " (" + data.batse.labels(2) + " < threshold)"] ...
          , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
        xlabel(data.batse.labels(1), "interpreter", "tex", "fontsize", fontSize);
        ylabel("dN / dT", "interpreter", "tex", "fontsize", fontSize);
        xlim([ax.xlower, ax.xupper]);
        xticks(10.^(linspace(ax.xlowTick, ax.xupTick, ax.xspan)));
        ylim([ax.ylower, ax.yupper]);
        yticks(10.^(linspace(ax.ylowTick, ax.yupTick, ax.yspan)));
    if saveNewImages
        export_fig("../out/Fig6_SynSam_dNdTvs" + xvar + ".png", "-m4 -transparent");
    end
hold off;

% Move figures to specific positions (in pixels) on the screen
movegui(fig1, [70, -1]);
movegui(fig2, [673, -1]);
movegui(fig3, [-70, -1]);
movegui(fig4, [70, 50]);
movegui(fig5, [673, 50]);
movegui(fig6, [-70, 50]);