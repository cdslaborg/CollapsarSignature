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
secondaryColor = '#F80'; %orange

% Options
histParam.T90.nbins = 45;
histParam.Epk.nbins = 45;
saveNewImages = false; % export figure on or off. Must have export_fig installed.

% T90 Data
tic
load("../out/synSamT90.mat"); % Load mat file
batse = importdata("..\in\BatseDur.xlsx");
icol.batse.T90 = 4; % column index of T90
data.batse.T90.data = batse.data(:, icol.batse.T90);
data.batse.T90.ndata = length(data.batse.T90.data);

histParam.T90.binEdges = getBinEdges(data.batse.T90.data, histParam.T90.nbins);
histParam.T90.binCenters = getBinCenters(histParam.T90.binEdges);
histParam.T90.binWidths = getBinWidths(histParam.T90.binEdges);
histParam.T90.countsVec = histcounts(data.batse.T90.data, histParam.T90.binEdges);
histParam.T90.logBinEdges = log(histParam.T90.binEdges);
histParam.T90.logBinCenters = log(histParam.T90.binCenters);

fo = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [44, -0.4, 1.7, 129, 3.5, 1.4]);
ft = fittype( ...
     'a / (c * sqrt(2 * pi)) * exp(-1 / 2 * ((x - b) / c)^2) + d / (g * sqrt(2 * pi)) * exp(-1 / 2 * ((x - f) / g)^2)' ...
     , 'options', fo); % fit a double Gaussian
data.batse.T90.fit = fit(histParam.T90.logBinCenters.', histParam.T90.countsVec.', ft);
toc

% Figure 1: Plot T90 histograms for BATSE and synSam
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', histParam.T90.binEdges, 'BinCounts', histParam.T90.countsVec);
    histogram('BinEdges', synSam.all.reduced.percentile50.binEdges, 'BinCounts', synSam.all.reduced.percentile50.countsVec, 'FaceColor', secondaryColor);
    plot(histParam.T90.binCenters, data.batse.T90.fit(histParam.T90.logBinCenters), 'color', 'black', 'lineWidth', lineWidth);
    plot(histParam.T90.binCenters, gaussian(data.batse.T90.fit.a, data.batse.T90.fit.b, data.batse.T90.fit.c, histParam.T90.logBinCenters) ...
        , 'color', 'blue', 'lineWidth', lineWidth);
    plot(histParam.T90.binCenters, gaussian(data.batse.T90.fit.d, data.batse.T90.fit.f, data.batse.T90.fit.g, histParam.T90.logBinCenters) ...
        , 'color', 'red', 'lineWidth', lineWidth);
%     xline( [min(data.batse.T90.data), min(synSam.all.binEdges), max(data.batse.T90.data) ...
%          , max(synSam.all.binEdges)], 'color', 'green');
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
    [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits([synSam.all.reduced.percentile50.binEdges, histParam.T90.binEdges], 'log');
    xlim([ax.lower, ax.upper]);
    xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    ylim([0, 200]);
    legend(["BATSE events, n = " + data.batse.T90.ndata + ", in " + histParam.T90.nbins + " bins" ...
          , "Synthetic events, n = " + synSam.all.reduced.percentile50.ndata + ", in " + length(synSam.all.reduced.percentile50.countsVec) + " bins" ...
          , "", "Fit: SGRBs", "Fit: LGRBs"], "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
hold off;

%{
x1 = normrnd(data.batse.f.b, data.batse.f.c, [1, 1000]);
x2 = normrnd(data.batse.f.f, data.batse.f.g, [1, 1000]);
expx1 = exp(x1);
expx2 = exp(x2);
binEdges = getBinEdges([expx1, expx2], nbins.T90);
binWidths = getBinWidths(binEdges);
countsVec = histcounts([expx1, expx2], binEdges);
countsVec1 = data.batse.f.a*histcounts(expx1, binEdges); ... data.batse.f.a*
countsVec2 = data.batse.f.d*histcounts(expx2, binEdges); ... data.batse.f.d*
countsVecS = countsVec1 + countsVec2;
area = sum(binWidths.*countsVecS);

figure("color", figureColor); hold on; box on;
    histogram('BinEdges', binEdges, 'BinCounts', countsVecS./area); ... ./280
    %plot(synSam.all.reduced.percentile50.binCenters, data.batse.f(logBinCenters), 'color', 'blue', 'lineWidth', lineWidth);
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
    [lower, upper, lowTick, upTick, span] = plotLimits(binEdges, 'log');
    xlim([lower, upper]);
    xticks(10.^(linspace(lowTick, upTick, span)));
    %ylim([0, 200]);
    title("Synthetic Data, n = " + sum(countsVec));
hold off;

figure("color", figureColor); hold on; box on;
    histogram('BinEdges', binEdges, 'BinCounts', countsVecS, 'normalization', 'probability');
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
    [lower, upper, lowTick, upTick, span] = plotLimits(binEdges, 'log');
    xlim([lower, upper]);
    xticks(10.^(linspace(lowTick, upTick, span)));
    title("Synthetic Data, n = " + sum(countsVec));
hold off;
%}







clear('synSam', 'batse', 'icol');

% Epk Data
tic
load("../out/synSamEpkT90.mat"); % Load mat file
synSam.both.all.detectedRand.ndata = length(synSam.both.all.detectedRand.Epk);
batse = importdata("..\in\Batse_orig.xlsx");
icol.batse.logEpk = 17;
data.batse.Epk.data = 10.^(batse.data.x1966GRBs(:, icol.batse.logEpk));
data.batse.Epk.ndata = length(data.batse.Epk.data);

%data.batse.Epk.binEdges = getBinEdges([data.batse.Epk.data.', synSam.both.all.detectedRand.Epk.'], histParam.Epk.nbins);
histParam.Epk.binEdges = getBinEdges(data.batse.Epk.data, histParam.Epk.nbins);
histParam.Epk.binCenters = getBinCenters(histParam.Epk.binEdges);
histParam.Epk.binWidths = getBinWidths(histParam.Epk.binEdges);
histParam.Epk.countsVec = histcounts(data.batse.Epk.data, histParam.Epk.binEdges);
synSam.both.all.detectedRand.countsVec = histcounts(synSam.both.all.detectedRand.Epk, histParam.Epk.binEdges);
histParam.Epk.logBinEdges = log(histParam.Epk.binEdges);
histParam.Epk.logBinCenters = log(histParam.Epk.binCenters);

fo = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [240, 5., 0.8, 44., 6.7, 0.5]);
ft = fittype( ...
     'a / (c * sqrt(2 * pi)) * exp(-1 / 2 * ((x - b) / c)^2) + d / (g * sqrt(2 * pi)) * exp(-1 / 2 * ((x - f) / g)^2)' ...
     , 'options', fo); % fit a double Gaussian
data.batse.Epk.fit = fit(histParam.Epk.logBinCenters.', histParam.Epk.countsVec.', ft);
toc

% Figure 2: Plot Epk histograms for BATSE and synSam
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', histParam.Epk.binEdges, 'BinCounts', histParam.Epk.countsVec);
    histogram('BinEdges', histParam.Epk.binEdges, 'BinCounts', synSam.both.all.detectedRand.countsVec, 'FaceColor', secondaryColor);
    plot(histParam.Epk.binCenters, data.batse.Epk.fit(histParam.Epk.logBinCenters), 'color', 'black', 'lineWidth', lineWidth);
    plot(histParam.Epk.binCenters, gaussian(data.batse.Epk.fit.a, data.batse.Epk.fit.b, data.batse.Epk.fit.c, histParam.Epk.logBinCenters) ...
        , 'color', 'red', 'lineWidth', lineWidth);
    plot(histParam.Epk.binCenters, gaussian(data.batse.Epk.fit.d, data.batse.Epk.fit.f, data.batse.Epk.fit.g, histParam.Epk.logBinCenters) ...
        , 'color', 'blue', 'lineWidth', lineWidth);
%     xline( [min(data.batse.Epk.data), min(synSam.both.all.detectedRand.Epk), max(data.batse.Epk.data) ...
%          , max(synSam.both.all.detectedRand.Epk)], 'color', 'green');
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
    xlabel("Epk [keV]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
    [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(histParam.Epk.binEdges, 'log');
    xlim([ax.lower, ax.upper]);
    xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    ylim([0, 200]);
    legend(["BATSE events, n = " + data.batse.Epk.ndata + ", in " + histParam.Epk.nbins + " bins" ...
          , "Synthetic events, n = " + synSam.both.all.detectedRand.ndata + ", in " + histParam.Epk.nbins + " bins" ...
          , "", "Fit: LGRBs", "Fit: SGRBs"], "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
hold off;

function gauss = gaussian(amp, mu, sigma, x)
    gauss = amp / (sigma * sqrt(2 * pi)) * exp(-1 / 2 * ((x - mu) / sigma).^2);
end