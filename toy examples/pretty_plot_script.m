addpath('~/Documents/MATLAB/export_fig-master')
addpath('~/Documents/MATLAB/prettyPlot')
clear all
load('1a.mat')
options.legendLoc = 'NorthEast';
options.logScale = 2;
options.colors = flip({'c','b','g',[1 0.5 0],'r'});
options.lineStyles =  {'-','-','-','-','-'};
options.markers = {'^','d','o','s','s'};
options.markerSpacing = [4000 1500
	4000 1000
	4000 500
	4000 800
	4000 200
	4000 500];
options.lineSize = 8;
% options.markerSpacing = markerSpacing;
% options.legendStr = names;
options.legend = Kappas;
options.ylabel = 'Duality Gap';
options.xlabel = 'Iteration ($\times 10^4$)';
options.legendFontSize=14;
options.labelLines = 0;
options.labelRotate = 1;
options.xlimits = [0 Kmax./10^4];
options.ylimits = [10^-6 5];
figure;
prettyPlot(fEvals,fVals,options);
export_fig('SP-FW','-pdf','-transparent');
clear all

load('1b.mat')
Kappas = {}
Kappas{end+1} = '$\nu = 0.09 \; \gamma =$ adaptive';
Kappas{end+1} = '$\nu = 0.09 \; \gamma =\frac{2}{2+k(t)}$';
Kappas{end+1} = '$\nu = 0.19 \; \gamma =$ adaptive';
Kappas{end+1} = '$\nu = 0.36 \; \gamma =$ adaptive';
Kappas{end+1} = '$\nu = 0.47 \; \gamma =$ adaptive';
sublin = fVals{1};
Esublin = fEvals{1};
fVals{1} = fVals{2};
fEvals{1} = fEvals{2};
fEvals{2} = Esublin;
fVals{2} = sublin;
options.legendLoc = 'NorthEast';
options.logScale = 2;
options.colors = flip({'c','b','g',[1 0.5 0],'r'});
options.lineStyles =  flip({'-','-','-','-','-'});
options.markers = {'^','d','o','s','s'};
options.markerSpacing = [4000 1500
	4000 1000
	4000 500
	4000 800
	4000 200
	4000 500];
options.lineSize = 8;
% options.markerSpacing = markerSpacing;
% options.legendStr = names;
options.legend = Kappas;
options.ylabel = 'Duality Gap';
options.xlabel = 'Iteration ($\times 10^4$)';
options.legendFontSize=14;
options.labelLines = 0;
options.labelRotate = 1;
options.xlimits = [0 Kmax./10^4];
options.ylimits = [10^-3 0.2];
figure;
prettyPlot(fEvals,fVals,options);
export_fig('SP-AFW','-pdf','-transparent');
clear all

load('1c.mat')
sublin = fVals{3};
Esublin = fEvals{3};
fVals{3} = fVals{4};
fEvals{3} = fEvals{4};
fEvals{4} = Esublin;
fVals{4} = sublin;
Kappas = {}
Kappas{end+1} = '$\nu = -4 \,\; \qquad  \gamma =$ heuristic';
Kappas{end+1} = '$\nu = -50 \qquad  \gamma =$ heuristic';
Kappas{end+1} = '$\nu = -5\cdot 10^4 \; \gamma =\frac{2}{2+k(t)}$';
Kappas{end+1} = '$\nu = -4\cdot 10^2  \; \gamma =$ heuristic';
Kappas{end+1} = '$\nu = -5\cdot 10^4  \; \gamma =$ heuristic';

options.legendLoc = 'NorthEast';
options.logScale = 2;
options.colors = flip({'c','b','g',[1 0.5 0],'r'});
options.lineStyles =  flip({'-','-','-','-','-'});
options.markers = {'^','d','o','s','s'};
options.markerSpacing = [2000 1500
	2000 200
	2000 500
	2000 800
	2000 400
	2000 500];
options.lineSize = 8;
% options.markerSpacing = markerSpacing;
% options.legendStr = names;
options.legend = flip(Kappas);
options.ylabel = 'Duality Gap'
options.xlabel = 'Iteration ($\times 10^4$)';


options.legendFontSize=14;
options.labelLines = 0;
options.labelRotate = 1;
options.xlimits = [0 Kmax./10^4];
options.ylimits = [10^-5 .1];
figure;
prettyPlot(flip(fEvals),flip(fVals),options);
export_fig('SP-AFW-heuristic','-pdf','-transparent');