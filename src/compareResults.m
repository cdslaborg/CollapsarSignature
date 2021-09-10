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
nbins.T90 = 45;
nbins.Epk = 45;
saveNewImages = false; % export figure on or off. Must have export_fig installed.

% T90 Data
tic
load("../out/synSamT90.mat"); % Load mat file
batse = importdata("..\in\BatseDur.xlsx");
icol.batse.T90 = 4; % column index of T90
data.batse.T90 = batse.data(:, icol.batse.T90);
data.batse.ndata = length(data.batse.T90);

binEdges = getBinEdges(data.batse.T90, nbins.T90);
binCenters = getBinCenters(binEdges);
binWidths = getBinCenters(binEdges);

data.batse.countsVec = histcounts(data.batse.T90, binEdges);
logBinEdges = log(binEdges);
logBinCenters = log(binCenters);
fo = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [44, -0.4, 1.7, 129, 3.5, 1.4]);
ft = fittype( ...
     'a / (c * sqrt(2 * pi)) * exp(-1 / 2 * ((x - b) / c)^2) + d / (g * sqrt(2 * pi)) * exp(-1 / 2 * ((x - f) / g)^2)' ...
     , 'options', fo); % fit a double Gaussian
data.batse.f = fit(logBinCenters.', data.batse.countsVec.', ft);
toc

% Figure 1: Plot T90 histograms for BATSE and synSam
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', binEdges, 'BinCounts', data.batse.countsVec);
    histogram('BinEdges', synSam.all.reduced.percentile50.binEdges, 'BinCounts', synSam.all.reduced.percentile50.countsVec, 'FaceColor', secondaryColor);
    plot(binCenters, data.batse.f(logBinCenters), 'color', 'blue', 'lineWidth', lineWidth);
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
    [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits([synSam.all.reduced.percentile50.binEdges, binEdges], 'log');
    xlim([ax.lower, ax.upper]);
    xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    ylim([0, 200]);
    legend(["BATSE events, n = " + data.batse.ndata + ", in " + nbins.T90 + " bins" ...
          , "Synthetic events, n = " + synSam.all.reduced.percentile50.ndata + ", in " + length(synSam.all.reduced.percentile50.countsVec) + " bins"] ...
          , "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
    %title("BATSE Data, n = " + data.batse.ndata);
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







clear('synSam', 'batse', 'icol', 'data');

% Epk Data
tic
load("../out/synSamEpkT90.mat"); % Load mat file
synSam.both.all.detectedRand.ndata = length(synSam.both.all.detectedRand.Epk);
batse = importdata("..\in\Batse_orig.xlsx");
icol.batse.logEpk = 17;
data.batse.Epk = 10.^(batse.data.x1966GRBs(:, icol.batse.logEpk));
data.batse.binEdges = getBinEdges([data.batse.Epk.', synSam.both.all.detectedRand.Epk.'], nbins.Epk);
data.batse.ndata = length(data.batse.Epk);
data.batse.countsVec = histcounts(data.batse.Epk, data.batse.binEdges);
synSam.both.all.detectedRand.countsVec = histcounts(synSam.both.all.detectedRand.Epk, data.batse.binEdges);
toc

% Figure 2: Plot Epk histograms for BATSE and synSam
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', data.batse.binEdges, 'BinCounts', data.batse.countsVec, "normalization", "probability");
    histogram('BinEdges', data.batse.binEdges, 'BinCounts', synSam.both.all.detectedRand.countsVec, "normalization", "probability", 'FaceColor', secondaryColor);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("Epk [keV]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
    [lower, upper, lowTick, upTick, span] = plotLimits(data.batse.binEdges, 'log');
    xlim([lower, upper]);
    xticks(10.^(linspace(lowTick, upTick, span)));
    ylim([5e-4, 4e-1]);
    legend(["BATSE events, n = " + data.batse.ndata + ", in " + nbins.Epk + " bins" ...
          , "Synthetic events, n = " + synSam.both.all.detectedRand.ndata + ", in " + nbins.Epk + " bins"] ...
          , "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
hold off;
