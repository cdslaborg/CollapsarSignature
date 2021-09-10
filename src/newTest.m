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
nbins.T50 = 43;
saveNewImages = false; % export figure on or off. Must have export_fig installed.

% T90 Data
tic
batse = importdata("..\in\BatseDur.xlsx");
batseSynT90 = importdata("..\in\batseSynthT90_process_1_sample.txt");
batseSynT50 = importdata("..\in\batseSynthT50_process_1_sample.txt");
icol.batseSat.T90 = 4; % column index of T90
icol.batseSat.T50 = 2;
icol.batseSyn.T90 = 2;
icol.batseSyn.T50 = 2;

data.batseSat.T90 = batse.data(:, icol.batseSat.T90);
data.batseSat.ndata = length(data.batseSat.T90);
data.batseSat.T50 = batse.data(:, icol.batseSat.T50);
data.batseSyn.T90 = exp(batseSynT90.data(:, icol.batseSyn.T90));
data.batseSyn.ndata = length(data.batseSyn.T90);
data.batseSyn.T50 = exp(batseSynT50.data(:, icol.batseSyn.T50));
clear('batse', 'batseSynT90', 'icol');

binEdges.T90 = getBinEdges([data.batseSat.T90; data.batseSyn.T90], nbins.T90);
binCenters.T90 = getBinCenters(binEdges.T90);
binWidths.T90 = getBinCenters(binEdges.T90);
binEdges.T50 = getBinEdges([data.batseSat.T50; data.batseSyn.T50], nbins.T50);
binCenters.T50 = getBinCenters(binEdges.T50);
binWidths.T50 = getBinCenters(binEdges.T50);

data.batseSat.T90countsVec = histcounts(data.batseSat.T90, binEdges.T90);
data.batseSat.T50countsVec = histcounts(data.batseSat.T50, binEdges.T50);
data.batseSyn.T90countsVec = histcounts(data.batseSyn.T90, binEdges.T90);
data.batseSyn.T50countsVec = histcounts(data.batseSyn.T50, binEdges.T50);

logBinEdges.T90 = log(binEdges.T90);
logBinCenters.T90 = log(binCenters.T90);
logBinEdges.T50 = log(binEdges.T50);
logBinCenters.T50 = log(binCenters.T50);

fo = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [193, 0.005, 1.7, 305, 3.6, 0.9]);
ft = fittype( ...
     'a / (c * sqrt(2 * pi)) * exp(-1 / 2 * ((x - b) / c)^2) + d / (g * sqrt(2 * pi)) * exp(-1 / 2 * ((x - f) / g)^2)' ...
     , 'options', fo); % fit a double Gaussian
data.batseSat.f.T90 = fit(logBinCenters.T90.', data.batseSat.T90countsVec.', ft);

fo = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [193, 0.005, 1.7, 305, 3.6, 0.9]);
ft = fittype( ...
     'a / (c * sqrt(2 * pi)) * exp(-1 / 2 * ((x - b) / c)^2) + d / (g * sqrt(2 * pi)) * exp(-1 / 2 * ((x - f) / g)^2)' ...
     , 'options', fo); % fit a double Gaussian
data.batseSat.f.T50 = fit(logBinCenters.T50.', data.batseSat.T50countsVec.', ft);
toc

% Figure 1: Plot T90 histograms for batseSat and batseSyn
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', binEdges.T90, 'BinCounts', data.batseSat.T90countsVec);
    histogram('BinEdges', binEdges.T90, 'BinCounts', data.batseSyn.T90countsVec, 'FaceColor', secondaryColor);
    plot(binCenters.T90, data.batseSat.f.T90(logBinCenters.T90), 'color', 'blue', 'lineWidth', lineWidth);
    xline([min(data.batseSat.T90), min(data.batseSyn.T90)], 'color', 'green');
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
    [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(binEdges.T90, 'log');
    xlim([ax.lower, ax.upper]);
    xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    ylim([0, 250]);
    legend(["BATSE events, n = " + data.batseSat.ndata + ", in " + nbins.T90 + " bins" ...
          , "New synthetic events, n = " + data.batseSyn.ndata + ", in " + nbins.T90 + " bins"] ...
          , "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
hold off;

% Figure 2: Plot T50 histograms for batseSat and batseSyn
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', binEdges.T50, 'BinCounts', data.batseSat.T50countsVec);
    histogram('BinEdges', binEdges.T50, 'BinCounts', data.batseSyn.T50countsVec, 'FaceColor', secondaryColor);
    plot(binCenters.T50, data.batseSat.f.T50(logBinCenters.T50), 'color', 'blue', 'lineWidth', lineWidth);
    xline([min(data.batseSat.T50), min(data.batseSyn.T50)], 'color', 'green');
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{50} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
    [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(binEdges.T50, 'log');
    xlim([ax.lower, ax.upper]);
    xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    ylim([0, 200]);
    legend(["BATSE events, n = " + data.batseSat.ndata + ", in " + nbins.T50 + " bins" ...
          , "New synthetic events, n = " + data.batseSyn.ndata + ", in " + nbins.T50 + " bins"] ...
          , "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
hold off;