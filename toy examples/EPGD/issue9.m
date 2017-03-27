% PEGD with SP in the interior
% with AFW projection Vs exact projection

% the library prettyplot needs to be loaded
%addpath('~/Documents/MATLAB/export_fig-master')
%addpath('~/Documents/MATLAB/prettyPlot')
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
c = 1
%===== Exact projection =====
alg = 0;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol,c);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD exact projection';
%===== 1 step for projection =====
alg = 1;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol,c);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD 1 projection';
%===== 3 step for projection =====
alg = 3;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol,c);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD 3 projection';
%===== 12 step for projection =====
alg = 12;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol,c);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD 12 projection';
%===== 24 step for projection =====
alg = 24;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol,c);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD 24 projection';
%===== \epsilon projection =====
alg = -1;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol,c);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD 1/t projection';
%===== \epsilon projection =====
alg = -1;
c = 1;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol,c);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD 5/t projection';
%===== \epsilon projection =====
alg = -1;
c = 5;
[G,niter] = EPGD(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol,c);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'EPGD 10/t projection';
%===== SPFW ======
adaptive = 4;
c = 10;
[G,niter] = SP_FW(Kmax,strong,M,cst,dim,a,b,alg,adaptive,tol);
fEvals{end+1} = (1:(niter-1));
fVals{end+1}=G(1:niter-1);
Mus{end+1} = 'SPFW';




% options.legendLoc = 'SouthWest';
options.logScale = 2;
options.colors = {   
	[0    0.4470    0.7410],
    [0.8500    0.3250    0.0980],
    [0.9290    0.6940    0.1250],
    [0.4940    0.1840    0.5560],
    [0.4660    0.6740    0.1880],
    [0.3010    0.7450    0.9330],
    [0.6350    0.0780    0.1840 ],
    [0    0.5000         0],
    [1.0000         0         0]
}
options.lineStyles =  {'--','--','-','-','-'};
options.markers ={'s','o','d','p','d'};
options.markerSpacing = [100 10
	1000 150
	1000 200
	1000 300
	1000 400
	1000 550
	1000 280
	1000 820
	1000 1000];
options.lineSize = 8;
% options.markerSpacing = markerSpacing;
% options.legendStr = names;
options.legend = Mus;
options.ylabel = 'Duality Gap';
options.xlabel = 'Oracle call';
options.legendFontSize=18;
options.labelLines = 0;
options.labelRotate = 1;
options.ylimits = [10^-5 10];
figure;
prettyPlot(fEvals,fVals,options);
export_fig('EPGD','-pdf','-transparent');

