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
histParam.T50.nbins = 43;
histParam.Epk.nbins = 45;
mergeBins = false; % whether to merge adjacent bins with <5 events
saveNewImages = false; % export figure on or off. Must have export_fig installed.

% T90 Data
tic
batse = importdata("..\in\BatseDur.xlsx");
batse2 = importdata("..\in\Batse_orig.xlsx");
batseSynT90 = importdata("..\in\batseSynthT90_process_1_sample.txt");
batseSynT50 = importdata("..\in\batseSynthT50_process_1_sample.txt");
batseSynEpk = importdata("..\in\batseSynthEpk_process_1_sample.txt");
icol.batseSat.T90 = 4; % column index of T90
icol.batseSat.T50 = 2;
icol.batseSat.logEpk = 17;
icol.batseSyn.T90 = 2;
icol.batseSyn.T50 = 2;
icol.batseSyn.Epk = 2;

data.batseSat.T90.data = batse.data(:, icol.batseSat.T90);
data.batseSat.T90.ndata = length(data.batseSat.T90.data);
data.batseSat.T50.data = batse.data(:, icol.batseSat.T50);
data.batseSat.T50.ndata = length(data.batseSat.T50.data);
data.batseSat.Epk.data = 10.^(batse2.data.x1966GRBs(:, icol.batseSat.logEpk));
data.batseSat.Epk.ndata = length(data.batseSat.Epk.data);

data.batseSyn.T90.data = exp(batseSynT90.data(:, icol.batseSyn.T90));
data.batseSyn.T90.ndata = length(data.batseSyn.T90.data);
data.batseSyn.T50.data = exp(batseSynT50.data(:, icol.batseSyn.T50));
data.batseSyn.T50.ndata = length(data.batseSyn.T50.data);
data.batseSyn.Epk.data = exp(batseSynEpk.data(:, icol.batseSyn.Epk));
data.batseSyn.Epk.ndata = length(data.batseSyn.Epk.data);
clear('batse', 'batse2', 'batseSynT90', 'batseSynT50', 'batseSynEpk', 'icol');

histParam.T90.binEdges = getBinEdges([data.batseSat.T90.data; data.batseSyn.T90.data], histParam.T90.nbins);
histParam.T90.binCenters = getBinCenters(histParam.T90.binEdges);
histParam.T90.binWidths = getBinWidths(histParam.T90.binEdges);
histParam.T50.binEdges = getBinEdges([data.batseSat.T50.data; data.batseSyn.T50.data], histParam.T50.nbins);
histParam.T50.binCenters = getBinCenters(histParam.T50.binEdges);
histParam.T50.binWidths = getBinWidths(histParam.T50.binEdges);
histParam.Epk.binEdges = getBinEdges([data.batseSat.Epk.data; data.batseSyn.Epk.data], histParam.Epk.nbins);
histParam.Epk.binCenters = getBinCenters(histParam.Epk.binEdges);
histParam.Epk.binWidths = getBinWidths(histParam.Epk.binEdges);

histParam.Sat.T90.countsVec = histcounts(data.batseSat.T90.data, histParam.T90.binEdges);
histParam.Sat.T50.countsVec = histcounts(data.batseSat.T50.data, histParam.T50.binEdges);
histParam.Sat.Epk.countsVec = histcounts(data.batseSat.Epk.data, histParam.Epk.binEdges);
histParam.Syn.T90.countsVec = histcounts(data.batseSyn.T90.data, histParam.T90.binEdges);
histParam.Syn.T50.countsVec = histcounts(data.batseSyn.T50.data, histParam.T50.binEdges);
histParam.Syn.Epk.countsVec = histcounts(data.batseSyn.Epk.data, histParam.Epk.binEdges);

histParam.T90.logBinEdges = log(histParam.T90.binEdges);
histParam.T90.logBinCenters = log(histParam.T90.binCenters);
histParam.T50.logBinEdges = log(histParam.T50.binEdges);
histParam.T50.logBinCenters = log(histParam.T50.binCenters);
histParam.Epk.logBinEdges = log(histParam.Epk.binEdges);
histParam.Epk.logBinCenters = log(histParam.Epk.binCenters);

% Fit Data to a double Gaussian
fun = 'a / (c * sqrt(2 * pi)) * exp(-1 / 2 * ((x - b) / c)^2) + d / (g * sqrt(2 * pi)) * exp(-1 / 2 * ((x - f) / g)^2)';
fo = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [193, 0.005, 1.7, 305, 3.6, 0.9]);
ft = fittype(fun, 'options', fo);
data.batseSat.T90.fit = fit(histParam.T90.logBinCenters.', histParam.Sat.T90.countsVec.', ft);
fo = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [193, 0.005, 1.7, 305, 3.6, 0.9]);
ft = fittype(fun, 'options', fo);
data.batseSat.T50.fit = fit(histParam.T50.logBinCenters.', histParam.Sat.T50.countsVec.', ft);
fo = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [240, 5., 0.8, 44., 6.7, 0.5]);
ft = fittype(fun, 'options', fo);
data.batseSat.Epk.fit = fit(histParam.Epk.logBinCenters.', histParam.Sat.Epk.countsVec.', ft);

% Convert the histogram counts to dN/dT (or dN/dE) by dividing each bin by its width
histParam.Sat.T90.dNdT = histParam.Sat.T90.countsVec ./ histParam.T90.binWidths;
histParam.Sat.T50.dNdT = histParam.Sat.T50.countsVec ./ histParam.T50.binWidths;
histParam.Sat.Epk.dNdE = histParam.Sat.Epk.countsVec ./ histParam.Epk.binWidths;
histParam.Syn.T90.dNdT = histParam.Syn.T90.countsVec ./ histParam.T90.binWidths;
histParam.Syn.T50.dNdT = histParam.Syn.T50.countsVec ./ histParam.T50.binWidths;
histParam.Syn.Epk.dNdE = histParam.Syn.Epk.countsVec ./ histParam.Epk.binWidths;
toc

% Figure 1: Plot T90 histograms for batseSat and batseSyn
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', histParam.T90.binEdges, 'BinCounts', histParam.Sat.T90.countsVec);
    histogram('BinEdges', histParam.T90.binEdges, 'BinCounts', histParam.Syn.T90.countsVec, 'FaceColor', secondaryColor);
    plot(histParam.T90.binCenters, data.batseSat.T90.fit(histParam.T90.logBinCenters), 'color', 'black', 'lineWidth', lineWidth);
    plot(histParam.T90.binCenters, gaussian(data.batseSat.T90.fit.a, data.batseSat.T90.fit.b, data.batseSat.T90.fit.c, histParam.T90.logBinCenters) ...
        , 'color', 'blue', 'lineWidth', lineWidth);
    plot(histParam.T90.binCenters, gaussian(data.batseSat.T90.fit.d, data.batseSat.T90.fit.f, data.batseSat.T90.fit.g, histParam.T90.logBinCenters) ...
        , 'color', 'red', 'lineWidth', lineWidth);
    %xline([min(data.batseSat.T90.data), min(data.batseSyn.T90.data)], 'color', 'green');
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
    [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(histParam.T90.binEdges, 'log');
    xlim([ax.lower, ax.upper]);
    xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    ylim([0, 230]);
    legend(["BATSE events, n = " + data.batseSat.T90.ndata + ", in " + histParam.T90.nbins + " bins" ...
          , "New synthetic events, n = " + data.batseSyn.T90.ndata + ", in " + histParam.T90.nbins + " bins" ...
          , "", "Fit: SGRBs", "Fit: LGRBs"], "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
hold off;

% Figure 2: Plot T50 histograms for batseSat and batseSyn
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', histParam.T50.binEdges, 'BinCounts', histParam.Sat.T50.countsVec);
    histogram('BinEdges', histParam.T50.binEdges, 'BinCounts', histParam.Syn.T50.countsVec, 'FaceColor', secondaryColor);
    plot(histParam.T50.binCenters, data.batseSat.T50.fit(histParam.T50.logBinCenters), 'color', 'black', 'lineWidth', lineWidth);
    plot(histParam.T50.binCenters, gaussian(data.batseSat.T50.fit.a, data.batseSat.T50.fit.b, data.batseSat.T50.fit.c, histParam.T50.logBinCenters) ...
        , 'color', 'blue', 'lineWidth', lineWidth);
    plot(histParam.T50.binCenters, gaussian(data.batseSat.T50.fit.d, data.batseSat.T50.fit.f, data.batseSat.T50.fit.g, histParam.T50.logBinCenters) ...
        , 'color', 'red', 'lineWidth', lineWidth);
    %xline([min(data.batseSat.T50.data), min(data.batseSyn.T50.data)], 'color', 'green');
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{50} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
    [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(histParam.T50.binEdges, 'log');
    xlim([ax.lower, ax.upper]);
    xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    ylim([0, 230]);
    legend(["BATSE events, n = " + data.batseSat.T90.ndata + ", in " + histParam.T50.nbins + " bins" ...
          , "New synthetic events, n = " + data.batseSyn.T90.ndata + ", in " + histParam.T50.nbins + " bins" ...
          , "", "Fit: SGRBs", "Fit: LGRBs"], "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
hold off;

% Figure 3: Plot Epk histograms for batseSat and batseSyn
figure("color", figureColor); hold on; box on;
    histogram('BinEdges', histParam.Epk.binEdges, 'BinCounts', histParam.Sat.Epk.countsVec);
    histogram('BinEdges', histParam.Epk.binEdges, 'BinCounts', histParam.Syn.Epk.countsVec, 'FaceColor', secondaryColor);
    plot(histParam.Epk.binCenters, data.batseSat.Epk.fit(histParam.Epk.logBinCenters), 'color', 'black', 'lineWidth', lineWidth);
    plot(histParam.Epk.binCenters, gaussian(data.batseSat.Epk.fit.a, data.batseSat.Epk.fit.b, data.batseSat.Epk.fit.c, histParam.Epk.logBinCenters) ...
        , 'color', 'red', 'lineWidth', lineWidth);
    plot(histParam.Epk.binCenters, gaussian(data.batseSat.Epk.fit.d, data.batseSat.Epk.fit.f, data.batseSat.Epk.fit.g, histParam.Epk.logBinCenters) ...
        , 'color', 'blue', 'lineWidth', lineWidth);
    set(gca, 'xscale', 'log', 'fontsize', fontSize-3);
    xlabel("E_{pk} [kEv]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("Counts", "interpreter", "tex", "fontsize", fontSize);
    [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(histParam.Epk.binEdges, 'log');
    xlim([ax.lower, ax.upper]);
    xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    ylim([0, 230]);
    legend(["BATSE events, n = " + data.batseSat.Epk.ndata + ", in " + histParam.Epk.nbins + " bins" ...
          , "New synthetic events, n = " + data.batseSyn.Epk.ndata + ", in " + histParam.Epk.nbins + " bins" ...
          , "", "Fit: LGRBs", "Fit: SGRBs"], "interpreter", "tex", "location", "northwest", "fontSize", fontSize-3);
hold off;

% Figure 4: Plot dNdT of T90
figure("color", figureColor); hold on; box on;
    if mergeBins
        [nbins1, ~, dNdT1] = binMergeStairPlot(histParam.Sat.T90.countsVec, histParam.T90.binEdges, defaultColor, lineWidth);
        [nbins2, ~, dNdT2] = binMergeStairPlot(histParam.Syn.T90.countsVec, histParam.T90.binEdges, secondaryColor, lineWidth);
        legend(["BATSE events, n = " + data.batseSat.T90.ndata + ", in " + nbins1 + " bins" ...
              , "New synthetic events, n = " + data.batseSyn.T90.ndata + ", in " + nbins2 + " bins"] ...
              , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits([dNdT1, dNdT2], 'log');
        ylim([ax.lower, ax.upper]);
        yticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    else
        hist2stairs(histParam.Sat.T90.dNdT, histParam.T90.binEdges, defaultColor, lineWidth);
        hist2stairs(histParam.Syn.T90.dNdT, histParam.T90.binEdges, secondaryColor, lineWidth);
        legend(["BATSE events, n = " + data.batseSat.Epk.ndata + ", in " + histParam.Epk.nbins + " bins" ...
              , "New synthetic events, n = " + data.batseSyn.Epk.ndata + ", in " + histParam.Epk.nbins + " bins"] ...
              , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits([histParam.Sat.T90.dNdT, histParam.Syn.T90.dNdT], 'log');
        ylim([ax.lower, ax.upper]);
        yticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT", "interpreter", "tex", "fontsize", fontSize);
    [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(histParam.T90.binEdges, 'log');
    xlim([ax.lower, ax.upper]);
    xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
hold off;

% Figure 5: Plot dNdT of T50
figure("color", figureColor); hold on; box on;
    if mergeBins
        [nbins1, ~, dNdT1] = binMergeStairPlot(histParam.Sat.T50.countsVec, histParam.T50.binEdges, defaultColor, lineWidth);
        [nbins2, ~, dNdT2] = binMergeStairPlot(histParam.Syn.T50.countsVec, histParam.T50.binEdges, secondaryColor, lineWidth);
        legend(["BATSE events, n = " + data.batseSat.T50.ndata + ", in " + nbins1 + " bins" ...
              , "New synthetic events, n = " + data.batseSyn.T50.ndata + ", in " + nbins2 + " bins"] ...
              , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits([dNdT1, dNdT2], 'log');
        ylim([ax.lower, ax.upper]);
        yticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    else
        hist2stairs(histParam.Sat.T50.dNdT, histParam.T50.binEdges, defaultColor, lineWidth);
        hist2stairs(histParam.Syn.T50.dNdT, histParam.T50.binEdges, secondaryColor, lineWidth);
        legend(["BATSE events, n = " + data.batseSat.T50.ndata + ", in " + histParam.T50.nbins + " bins" ...
              , "New synthetic events, n = " + data.batseSyn.T50.ndata + ", in " + histParam.T50.nbins + " bins"] ...
              , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
        [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits([histParam.Sat.T50.dNdT, histParam.Syn.T50.dNdT], 'log');
        ylim([ax.lower, ax.upper]);
        yticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
    end
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{50} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT", "interpreter", "tex", "fontsize", fontSize);
    [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(histParam.T50.binEdges, 'log');
    xlim([ax.lower, ax.upper]);
    xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
hold off;

% Figure 6: Plot dNdE of Epk
% figure("color", figureColor); hold on; box on;
%     if mergeBins
%         [nbins1, ~, dNdE1] = binMergeStairPlot(histParam.Sat.Epk.countsVec, histParam.Epk.binEdges, defaultColor, lineWidth);
%         [nbins2, ~, dNdE2] = binMergeStairPlot(histParam.Syn.Epk.countsVec, histParam.Epk.binEdges, secondaryColor, lineWidth);
%         legend(["BATSE events, n = " + data.batseSat.Epk.ndata + ", in " + nbins1 + " bins" ...
%               , "New synthetic events, n = " + data.batseSyn.Epk.ndata + ", in " + nbins2 + " bins"] ...
%               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
%         [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits([dNdE1, dNdE2], 'log');
%         ylim([ax.lower, ax.upper]);
%         yticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
%     else
%         hist2stairs(histParam.Sat.Epk.dNdE, histParam.Epk.binEdges, defaultColor, lineWidth);
%         hist2stairs(histParam.Syn.Epk.dNdE, histParam.Epk.binEdges, secondaryColor, lineWidth);
%         legend(["BATSE events, n = " + data.batseSat.Epk.ndata + ", in " + histParam.Epk.nbins + " bins" ...
%               , "New synthetic events, n = " + data.batseSyn.Epk.ndata + ", in " + histParam.Epk.nbins + " bins"] ...
%               , "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
%         [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits([histParam.Sat.Epk.dNdE, histParam.Syn.Epk.dNdE], 'log');
%         ylim([ax.lower, ax.upper]);
%         yticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
%     end
%     set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
%     xlabel("E_{pk} [kEv]", "interpreter", "tex", "fontsize", fontSize);
%     ylabel("dN / dE", "interpreter", "tex", "fontsize", fontSize);
%     [ax.lower, ax.upper, ax.lowTick, ax.upTick, ax.span] = plotLimits(histParam.Epk.binEdges, 'log');
%     xlim([ax.lower, ax.upper]);
%     xticks(10.^(linspace(ax.lowTick, ax.upTick, ax.span)));
% hold off;


function gauss = gaussian(amp, mu, sigma, x)
    gauss = amp / (sigma * sqrt(2 * pi)) * exp(-1 / 2 * ((x - mu) / sigma).^2);
end