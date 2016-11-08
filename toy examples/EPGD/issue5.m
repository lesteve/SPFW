% PEGD with SP in the interior
% with AFW projection Vs exact projection

% the library prettyplot needs to be loaded
addpath('~/Documents/MATLAB/export_fig-master')
addpath('~/Documents/MATLAB/prettyPlot')
% seed for reproductibility
rng(2)

% ==== INITIALIZATION ====
fEvals={};
fVals= {};
Mus = {};
clear options
strong = 1; %FOR SP-AFW
cst = .1;
dim = 30;
a = [rand(dim,1)./2 + .25]; 
b = [rand(dim,1)./2 + .25]; 
M = cst.* (2.*rand(dim,dim)-1); 
Kmax = 4000;
tol = 10^(-5);

%===== constant step =====
adaptive = 0;
%===== Exact projection =====
alg = 0;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD exact projection';
%===== 1 step for projection =====
alg = 1;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD 1 projection';
%===== 3 step for projection =====
alg = 3;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD 3 projection';
%===== \epsilon projection =====
alg = -1;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD gap projection';
%===== SPFW ======
adaptive = 4;
[G,niter] = SP_FW(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'SPFW';




% options.legendLoc = 'SouthWest';
options.logScale = 2;
options.colors = {'c','b','g',[1 0.5 0],'r'};
options.lineStyles =  {'--','--','-','-','-'};
options.markers ={'s','o','d','p','d'};
options.markerSpacing = [400 100
	400 1000
	400 50
	400 80
	400 20
	400 50];
options.lineSize = 8;
% options.markerSpacing = markerSpacing;
% options.legendStr = names;
options.legend = Mus;
options.ylabel = 'Duality Gap';
options.xlabel = 'Iteration';
options.legendFontSize=18;
options.labelLines = 0;
options.labelRotate = 1;
figure;
prettyPlot(fEvals,fVals,options);