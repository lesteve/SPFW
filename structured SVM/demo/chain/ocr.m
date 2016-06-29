% Applies the structured SVM to the OCR dataset by Ben Taskar. The structured
% model considered here is the standard chain graph, with the pixel values of
% the digit as unary features and a transition matrix of size num_states^2 as
% a pairwise potential. Additionally, we include a unary bias term for the first
% and last symbol in the sequence.

addpath(genpath('../../solvers/'));
addpath('helpers');

% We support two different settings for the dataset (ocr: only one fold in
% training set, ocr2: all but one fold in training set
% -- ocr2 is the one that we have used in our experiments in the 
% ICML 2013 paper)
data_name = 'ocr';
[patterns_train, labels_train, patterns_test, labels_test] = loadOCRData(data_name, '../../data/');

%% == run one of the solvers on the problem

% create problem structure:
subset_size = 100;
param = [];
param.patterns = patterns_train(1:subset_size);
param.labels = labels_train(1:subset_size);
param.lossFn = @chain_loss;
param.oracleFn = @chain_oracle;
param.featureFn = @chain_featuremap;

% options structure:
options = [];
options.alpha = 2;
options.lambda = 0.0;
options.gap_threshold = 0.00001; % duality gap stopping criterion
options.num_passes = 2000 ; % max number of passes through data
options.do_line_search = 0;
options.sample = 'uniform'; 
options.debug_iter = 100;
beta = 0.01;

% Compute the approximal optimal value
options.solution = 0;
options.debug = 1;
options.average = 0;
options.beta = BB(1);

param.stepsize =1;
[model, progress] = solverSP_BCFW(param, options);
w_star = model.w;
% save('w_star')
options.solution = 1;

options.w_star = w_star;


options.num_passes = 100 ; % max number of passes through data

options.beta = beta;
[model1,progress1] = solverSubgradient(param,options);
[model2,progress2] = solverStoSubgradient(param,options);
param.stepsize=0.1;
[model3,progress3] = solverStoSubgradient(param,options);
avg_loss = 0;
% for i=1:subset_size
%   ypredict = chain_oracle(param, model3, patterns_train{i}); % standard prediction as don't give label as input
%   avg_loss = avg_loss + chain_loss(param, labels_train{i}, ypredict);
% end
% avg_loss = avg_loss / subset_size;
% fprintf('average loss on the training set: %f.\n', avg_loss);

% % loss on test set
% avg_loss = 0;
% for i=1:numel(patterns_test)
%     ypredict = chain_oracle(param, model3, patterns_test{i});
%     avg_loss = avg_loss + chain_loss(param, labels_test{i}, ypredict);
% end
% avg_loss = avg_loss / numel(patterns_test);
% fprintf('average loss on the test set: %f.\n', avg_loss);
[model4,progress4] = solverSP_FW(param,options);
options.alpha = 1;
[model5,progress5] = solverSP_FW(param,options);
options.average =1;
[model6,progress6] = solverSP_BCFW(param,options);
loglog(progress1.eff_pass, progress1.primal,'b'); % primal
hold on
loglog(progress2.eff_pass, progress2.primal,'c')
loglog(progress3.eff_pass, progress3.primal,'g'); % primal
loglog(progress4.eff_pass, progress4.primal,'r'); % primal
loglog(progress5.eff_pass, progress5.primal,'c'); % primal
loglog(progress6.eff_pass, progress6.primal,'y'); % primal
legend('Subgradient','SSG','SSG2','SPFW','SPFW2','BCSWFW')
hold off
