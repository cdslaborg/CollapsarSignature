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
numBinsBatse = 50; % number of bins to bin BATSE data in
numBinsSwift = 50; % number of bins to bin Swift data in
numBinsFermi = 50; % number of bins to bin Fermi data in


% Begin Data Processing
batse = importdata("..\in\BATSE.xlsx"); % Complete dataset
swift = importdata("..\in\Swift.txt"); % Data aquired 2021-07-11
swift.data = rmmissing(swift.data); % remove NaN values
fermi = importdata("..\in\Fermi.txt"); % Data aquired 2021-07-11
fermi.data = rmmissing(fermi.data); % remove NaN values

icol = struct();
icol.batse.T90 = 14; % column index of T90
icol.swift.T90 = 1; % column index of T90
icol.fermi.T90 = 1; % column index of T90

data = struct(); 
data.batse.T90 = batse.data.x1966GRBs(:,icol.batse.T90);
data.batse.binEdges = getBinEdges(data.batse.T90,numBinsBatse);
data.batse.binCenters = getBinCenters(data.batse.binEdges);
data.swift.T90 = swift.data(:,icol.swift.T90);
data.swift.binEdges = getBinEdges(data.swift.T90,numBinsSwift);
data.swift.binCenters = getBinCenters(data.swift.binEdges);
data.fermi.T90 = fermi.data(:,icol.fermi.T90);
data.fermi.binEdges = getBinEdges(data.fermi.T90,numBinsFermi);
data.fermi.binCenters = getBinCenters(data.fermi.binEdges);

data.batse.ndata = length(data.batse.T90);
data.batse.countsVec = histcounts(data.batse.T90,data.batse.binEdges);
data.swift.ndata = length(data.swift.T90);
data.swift.countsVec = histcounts(data.swift.T90,data.swift.binEdges);
data.fermi.ndata = length(data.fermi.T90);
data.fermi.countsVec = histcounts(data.fermi.T90,data.fermi.binEdges);
clear('batse','swift','fermi','icol');

% Convert the dependent variable dN/dlog(T) to dN/dT by dividing it by T90
data.batse.dNdT = data.batse.countsVec ./ data.batse.binCenters;
data.swift.dNdT = data.swift.countsVec ./ data.swift.binCenters;
data.fermi.dNdT = data.fermi.countsVec ./ data.fermi.binCenters;

% End Data Processing


% Figure 1: Plot dN/dT of BATSE dataset
figure("color", figureColor); hold on; box on;
    hist2stairs(data.batse.dNdT, data.batse.binEdges, defaultColor, lineWidth);
    [nbins, binEdges, dNdT] = binMergeStairPlot(data.batse.countsVec, data.batse.binEdges, secondaryColor, lineWidth);
    %xline([2e-2, 2e2], 'color', 'blue', 'lineWidth', lineWidth);
    [lower, upper] = plotLimits(binEdges, 'log');
    xlim([lower, upper]);
    [lower, upper] = plotLimits([data.batse.dNdT, dNdT], 'log');
    ylim([lower, upper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
    legend([data.batse.ndata + " BATSE observed GRBs in " + numBinsBatse + " bins", "Adjacent bins with <5 events merged," + newline ...
           + "resulting in " + nbins + " bins"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    if saveNewImages
        export_fig(data.output.path + "/" + "BatseDNDT.png", "-m4 -transparent");
    end
hold off;

% Figure 2: Plot dN/dT of Swift dataset
figure("color", figureColor); hold on; box on;
    hist2stairs(data.swift.dNdT, data.swift.binEdges, defaultColor, lineWidth);
    [nbins, binEdges, dNdT] = binMergeStairPlot(data.swift.countsVec, data.swift.binEdges, secondaryColor, lineWidth);
    %xline([4e-2, 2e2], 'color', 'magenta', 'lineWidth', lineWidth);
    [lower, upper] = plotLimits(binEdges, 'log');
    xlim([lower, upper]);
    [lower, upper] = plotLimits([data.swift.dNdT, dNdT], 'log');
    ylim([lower, upper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
    legend([data.swift.ndata + " Swift observed GRBs in " + numBinsSwift + " bins", "Adjacent bins with <5 events merged," + newline ...
           + "resulting in " + nbins + " bins"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    if saveNewImages
        export_fig(data.output.path + "/" + "SwiftDNDT.png", "-m4 -transparent");
    end
hold off;

% Figure 3: Plot dN/dT of Fermi dataset
figure("color", figureColor); hold on; box on;
    hist2stairs(data.fermi.dNdT, data.fermi.binEdges, defaultColor, lineWidth);
    [nbins, binEdges, dNdT] = binMergeStairPlot(data.fermi.countsVec, data.fermi.binEdges, secondaryColor, lineWidth);
    %xline([1e-1, 2e2], 'color', 'green', 'lineWidth', lineWidth);
    [lower, upper] = plotLimits(binEdges, 'log');
    xlim([lower, upper]);
    [lower, upper] = plotLimits([data.fermi.dNdT, dNdT], 'log');
    ylim([lower, upper]);
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fontSize-3);
    xlabel("T_{90} [s]", "interpreter", "tex", "fontsize", fontSize);
    ylabel("dN / dT_{90} [s^{-1}]", "interpreter", "tex", "fontsize", fontSize);
    legend([data.fermi.ndata + " Fermi observed GRBs in " + numBinsFermi + " bins", "Adjacent bins with <5 events merged," + newline ...
           + "resulting in " + nbins + " bins"], "interpreter", "tex", "location", "southwest", "fontSize", fontSize-3);
    if saveNewImages
        export_fig(data.output.path + "/" + "FermiDNDT.png", "-m4 -transparent");
    end
hold off;
