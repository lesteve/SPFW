% the library prettyplot needs to be loaded
addpath('~/Documents/MATLAB/export_fig-master')
addpath('~/Documents/MATLAB/prettyPlot')
% seed for reproductibility
rng(2)

% ==== INITIALIZATION ====
fEvals={};
fVals= {};
taus = {};
clear options
strongs = [10;2;0.5;0.1]; %FOR SP-AFW
cst = .1;
dim = 30;
a_away = [-(sign(2*rand(dim/2,1)-1)-1)*.5; rand(dim/2,1)]; 
b_away = [-(sign(2*rand(dim/2,1)-1)-1)*.5; rand(dim/2,1)]; 
a_int =  rand(dim,1); 
b_int =  rand(dim,1);
M = cst.* (2.*rand(dim,dim)-1); 
adaptive = 3;
Kmax = 10000;
for i = 1:4
	adaptive = 3;
	strong = strongs(i);
	away = 1;
	[G,niter,tau] = SP_FW(Kmax,strong,M,cst,dim,a_away,b_away,away,adaptive);
	fEvals{end+1} = 1:(niter-1);
	fVals{end+1}=G(1:niter-1);
	taus{end+1} = strcat('SP-AFW');
	away = 0;
	[G,niter,tau] = SP_FW(Kmax,strong,M,cst,dim,a_int,b_int,away,adaptive);
	fEvals{end+1} = 1:(niter-1);
	fVals{end+1}=G(1:niter-1);
	taus{end+1} = strcat('SP-FW');
end
fEvals{end+1} = 1:(niter-1);
fVals{end+1}=1./(1:niter-1);
fEvals{end+1} = 1:(niter-1);
fVals{end+1}=1./sqrt((1:niter-1));
% options.legendLoc = 'SouthWest';
options.logScale = 3;
options.colors = {'c','c','b','b','g','g',[1 0.5 0],[1 0.5 0],'r'};
options.lineStyles =  {'--','-','--','-','--','-','--','-','-'};
% options.markers ={'s','s','o','o','d','d','p','p','d'};
% options.markerSpacing = [2000 150
% 	2000 20
% 	2000 50
% 	2000 80
% 	2000 50
% 	2000 80
% 	2000 40
% 	2000 50];
options.lineSize = 8;
% options.markerSpacing = markerSpacing;
% options.legendStr = names;
options.legend =taus;
options.ylabel = 'Duality Gap';
options.xlabel = 'Iteration';
options.legendFontSize=18;
options.labelLines = 0;
options.labelRotate = 1;
options.xlimits = [0 Kmax];
% options.ylimits = [10^-3 1];
figure;
prettyPlot(fEvals,fVals,options);