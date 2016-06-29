% the library prettyplot needs to be loaded

fEvals={};
fVals= {};
alphas = {};
clear options
% strongs = [37.6 ; 40; 50;100]
% strongs = [8 ; 4;1;0.1]; % for alpha neg
strongs = [8.3;10;15;80]; %FOR SP-AFW
cst = .1;
dim = 30;
%set solution at the boundary
a = 2*ones(dim,1)-1; b=2*ones(dim,1)-1; 

% strongs = [1.15 ;1.5; 3;10]; %For int-SPFW
% cst = .04;
% dim = 5 
adaptive = 3; %Cf SP_FW code
Kmax = 2000;
strong = strongs(1);
SP_FW
fEvals{end+1} = 1:(iter-1);
fVals{end+1}=G(1:iter-1);
alphas{end+1} = strcat('\it\alpha-1/2=',num2str(alpha,2),'\rm \gamma = ',' 2/(2+k)');
for i = 1:4
	adaptive = 1
	strong = strongs(i);
	SP_FW
	fEvals{end+1} = 1:(iter-1);
	fVals{end+1}=G(1:iter-1)
	i
	if i == 1
		alphas{end+1} = strcat('\it\alpha=',num2str(alpha,2),'\rm \gamma adaptive');
	else
		alphas{end+1} = strcat('\alpha=',num2str(alpha,2),'\rm   \gamma adaptive');
	end
end

options.legendLoc = 'SouthWest';
options.logScale = 2;
options.colors = {'c','b','g',[1 0.5 0],'r'};
options.lineStyles =  {'--','--','-','-','-'};
options.markers ={'s','o','d','p','d'};
options.markerSpacing = [200 150
	200 100
	100 50
	100 80
	100 20
	100 50];
options.lineSize = 8;
% options.markerSpacing = markerSpacing;
% options.legendStr = names;
options.legend =alphas;
options.ylabel = 'Duality Gap';
options.xlabel = 'Iteration';
options.legendFontSize=18;
options.labelLines = 0;
options.labelRotate = 1;
options.xlimits = [100 2000];
options.ylimits = [10^-5 1];
figure;
prettyPlot(fEvals,fVals,options);
% legend('SPFW','BCSWFW','Subgradient','SSG','SSG')
% print -dpdf prettyPlot4.pdf