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

% Options
saveNewImages = false; % export figure on or off. Must have export_fig installed.
numBinsBatse = 43; % number of bins to bin BATSE data in
numBinsSwift = 40; % number of bins to bin Swift data in
numBinsFermi = 40; % number of bins to bin Fermi data in


% Begin Data Processing (data aquired 2021-07-11)
tic
batse = importdata("..\in\BatseDur.xlsx");
swift = importdata("..\in\Swift.txt");
swift.data = rmmissing(swift.data); % remove NaN values
fermi = importdata("..\in\Fermi.txt");
fermi.data = rmmissing(fermi.data); % remove NaN values

icol = struct();
icol.batse.T50 = 2; % column index of T50
icol.batse.T90 = 4; % column index of T90
icol.swift.T90 = 1;
icol.fermi.T90 = 1;

data = struct(); 
data.batse.T90 = batse.data(:, icol.batse.T90);
data.batse.binEdges = getBinEdges(data.batse.T90, numBinsBatse);
data.batse.binWidths = getBinWidths(data.batse.binEdges);
data.batse.T50.data = batse.data(:, icol.batse.T50);
data.batse.T50.binEdges = getBinEdges(data.batse.T50.data, numBinsBatse);
data.batse.T50.binWidths = getBinWidths(data.batse.T50.binEdges);
data.swift.T90 = swift.data(:, icol.swift.T90);
data.swift.binEdges = getBinEdges(data.swift.T90, numBinsSwift);
data.swift.binWidths = getBinWidths(data.swift.binEdges);
data.fermi.T90 = fermi.data(:, icol.fermi.T90);
data.fermi.binEdges = getBinEdges(data.fermi.T90, numBinsFermi);
data.fermi.binWidths = getBinWidths(data.fermi.binEdges);

data.batse.ndata = length(data.batse.T90);
data.batse.countsVec = histcounts(data.batse.T90, data.batse.binEdges);
data.batse.T50.countsVec = histcounts(data.batse.T50.data, data.batse.T50.binEdges);
data.swift.ndata = length(data.swift.T90);
data.swift.countsVec = histcounts(data.swift.T90, data.swift.binEdges);
data.fermi.ndata = length(data.fermi.T90);
data.fermi.countsVec = histcounts(data.fermi.T90, data.fermi.binEdges);
clear('batse','swift','fermi','icol');

% Convert the dependent variable dN/dlog(T) to dN/dT by dividing it by T90
data.batse.dNdT = data.batse.countsVec ./ data.batse.binWidths;
data.batse.T50.dNdT = data.batse.T50.countsVec ./ data.batse.T50.binWidths;
data.swift.dNdT = data.swift.countsVec ./ data.swift.binWidths;
data.fermi.dNdT = data.fermi.countsVec ./ data.fermi.binWidths;
toc
% End Data Processing

% Figure 1: Plot dN/dT of BATSE T90 dataset
figure("color", figureColor); hold on; box on;
    hist2stairs(data.batse.dNdT, data.batse.binEdges, defaultColor, lineWidth);
    [nbins, binEdges, dNdT] = binMergeStairPlot(data.batse.countsVec, data.batse.binEdges, secondaryColor, lineWidth);
    %xline([2e-2, 2e2], 'color', 'blue', 'lineWidth', lineWidth);
    [lower, upper, lowTick, upTick, span] = plotLimits(binEdges, 'log');
    xlim([lower, upper]);
    xticks(10.^(linspace(lowTick, upTick, span)));
    [lower, upper, lowTick, upTick, span] = plotLimits([data.batse.dNdT, dNdT], 'log');
    ylim([lower, upper]);
    yticks(10.^(linspace(lowTick, upTick, span)));
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
    legend([data.batse.ndata + " BATSE observed GRBs in " + numBinsBatse + " bins", "Adjacent bins with <5 events merged," + newline ...
           + "resulting in " + nbins + " bins"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    if saveNewImages
        export_fig(data.output.path + "/" + "BatseDNDT.png", "-m4 -transparent");
    end
hold off;

% Figure 2: Plot dN/dT of BATSE T50 dataset
figure("color", figureColor); hold on; box on;
    hist2stairs(data.batse.T50.dNdT, data.batse.T50.binEdges, defaultColor, lineWidth);
    [nbins, binEdges, dNdT] = binMergeStairPlot(data.batse.T50.countsVec, data.batse.T50.binEdges, secondaryColor, lineWidth);
    [lower, upper, lowTick, upTick, span] = plotLimits(binEdges, 'log');
    xlim([lower, upper]);
    xticks(10.^(linspace(lowTick, upTick, span)));
    [lower, upper, lowTick, upTick, span] = plotLimits([data.batse.T50.dNdT, dNdT], 'log');
    ylim([lower, upper]);
    yticks(10.^(linspace(lowTick, upTick, span)));
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{50} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{50} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
    legend([data.batse.ndata + " BATSE observed GRBs in " + numBinsBatse + " bins", "Adjacent bins with <5 events merged," + newline ...
           + "resulting in " + nbins + " bins"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    if saveNewImages
        export_fig(data.output.path + "/" + "BatseT50DNDT.png", "-m4 -transparent");
    end
hold off;

% Figure 3: Plot dN/dT of Swift dataset
figure("color", figureColor); hold on; box on;
    hist2stairs(data.swift.dNdT, data.swift.binEdges, defaultColor, lineWidth);
    [nbins, binEdges, dNdT] = binMergeStairPlot(data.swift.countsVec, data.swift.binEdges, secondaryColor, lineWidth);
    %xline([4e-2, 2e2], 'color', 'magenta', 'lineWidth', lineWidth);
    [lower, upper, lowTick, upTick, span] = plotLimits(binEdges, 'log');
    xlim([lower, upper]);
    xticks(10.^(linspace(lowTick, upTick, span)));
    [lower, upper, lowTick, upTick, span] = plotLimits([data.swift.dNdT, dNdT], 'log');
    ylim([lower, upper]);
    yticks(10.^(linspace(lowTick, upTick, span)));
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
    legend([data.swift.ndata + " Swift observed GRBs in " + numBinsSwift + " bins", "Adjacent bins with <5 events merged," + newline ...
           + "resulting in " + nbins + " bins"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    if saveNewImages
        export_fig(data.output.path + "/" + "SwiftDNDT.png", "-m4 -transparent");
    end
hold off;

% Figure 4: Plot dN/dT of Fermi dataset
figure("color", figureColor); hold on; box on;
    hist2stairs(data.fermi.dNdT, data.fermi.binEdges, defaultColor, lineWidth);
    [nbins, binEdges, dNdT] = binMergeStairPlot(data.fermi.countsVec, data.fermi.binEdges, secondaryColor, lineWidth);
    %xline([1e-1, 2e2], 'color', 'green', 'lineWidth', lineWidth);
    [lower, upper, lowTick, upTick, span] = plotLimits(binEdges, 'log');
    xlim([lower, upper]);
    xticks(10.^(linspace(lowTick, upTick, span)));
    [lower, upper, lowTick, upTick, span] = plotLimits([data.fermi.dNdT, dNdT], 'log');
    ylim([lower, upper]);
    yticks(10.^(linspace(lowTick, upTick, span)));
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
    legend([data.fermi.ndata + " Fermi observed GRBs in " + numBinsFermi + " bins", "Adjacent bins with <5 events merged," + newline ...
           + "resulting in " + nbins + " bins"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    if saveNewImages
        export_fig(data.output.path + "/" + "FermiDNDT.png", "-m4 -transparent");
    end
hold off;