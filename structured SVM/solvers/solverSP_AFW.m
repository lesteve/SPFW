function [model, progress] = solverSPFW(param, options)
% [model, progress] = solverFW(param, options)
%
% Solves the structured support vector machine (SVM) using the slower batch
% Frank-Wolfe algorithm, see (Lacoste-Julien, Jaggi, Schmidt, Pletscher; ICML
% 2013) for more details.
% This is Algorithm 2 in the paper, and the code here follows a similar
% notation. Each step of the method calls the decoding oracle (i.e. the 
% loss-augmented predicion) for all points.
%
% The structured SVM has the form
% min_{w} 0.5*\lambda*||w||^2+ 1/n*\sum_{i=1}^n H_i(w) [see (3) in paper]
%   where H_i(w) is the structured hinge loss on the i^th example:
%         H_i(w) = max_{y in Y} L(y_i,y) - <w, \psi_i(y)> [(2) in paper]
%
% We use a calling interface very similar to version 1.1 of svm-struct-matlab 
% developped by Andrea Vedaldi (see vedaldi/code/svm-struct-matlab.html).
% svm-struct-matlab is a Matlab wrapper interface to the widely used
% SVM^struct code by Thorsten Joachims (http://svmlight.joachims.org/svm_struct.html) 
% which implements a cutting plane algorithm. 
% 
% If your code was using:
%    model = svm_struct_learn(command_args,param)
%
% You can replace it with:
%    model = solverFW(param, options)
% 
% with the same param structure and letting options.lambda = 1/C
% to solve the *same* optimization problem. [Make sure your code only 
% uses model.w as we currently don't return the dual variables, etc. in the
% model structure unlike svm-struct-learn].
% 
% Inputs:
%   param: a structure describing the problem with the following fields:
%
%     patterns  -- patterns (x_i)
%         A cell array of patterns (x_i). The entries can have any
%         nature (they can just be indexes of the actual data for
%         example).
%     
%     labels    -- labels (y_i)
%         A cell array of labels (y_i). The entries can have any nature.
%
%     lossFn    -- loss function callback
%         A handle to the loss function L(ytruth, ypredict) defined for 
%         your problem. This function should have a signature of the form:
%           scalar_output = loss(param, ytruth, ypredict) 
%         It will be given an input ytruth, a ground truth label;
%         ypredict, a prediction label; and param, the same structure 
%         passed to solverFW.
% 
%     oracleFn  -- loss-augmented decoding callback
%         [Can also be called constraintFn for backward compatibility with
%          code using svm_struct_learn.]
%         A handle to the 'maximization oracle' (equation (2) in paper) 
%         which solves the loss-augmented decoding problem. This function
%         should have a signature of the form:
%           ypredict = decode(param, model, x, y)
%         where x is an input pattern, y is its ground truth label,
%         param is the input param structure to solverFW and model is the
%         current model structure (the main field is model.w which contains
%         the parameter vector).
%
%     featureFn  feature map callback
%         A handle to the feature map function \phi(x,y). This function
%         should have a signature of the form:
%           phi_vector = feature(param, x, y)
%         where x is an input pattern, y is an input label, and param 
%         is the usual input param structure. The output should be a vector 
%         of *fixed* dimension d which is the same
%         across all calls to the function. The parameter vector w will
%         have the same dimension as this feature map. In our current
%         implementation, w is sparse if phi_vector is sparse.
% 
%  options:    (an optional) structure with some of the following fields to
%              customize the behavior of the optimization algorithm:
% 
%   lambda      The regularization constant (default: 1/n).
%
%   gap_threshold **STOPPING CRITERION**
%               Stop the algorithm once the duality gap falls below this
%               threshold. Note that the default of 0.1 is equivalent
%               to the criterion used in svm_struct_learn Matlab wrapper.
%               (default: 0.1).
%
%   num_passes  Maximum number of passes through the data before the 
%               algorithm stops (default: 200)
%
%   debug       Boolean flag whether to track the primal objective, dual
%               objective, and training error (makes the code about 3x
%               slower given the extra two passes through data).
%               (default: 0)
%   do_linesearch
%               Boolean flag whether to use line-search. (default: 1)
%   time_budget Number of minutes after which the algorithm should terminate.
%               Useful if the solver is run on a cluster with some runtime
%               limits. (default: inf)
%   test_data   Struct with two fields: patterns and labels, which are cell
%               arrays of the same form as the training data. If provided the
%               logging will also evaluate the test error.
%               (default: [])
%
% Outputs:
%   model       model.w contains the parameters;
%               model.ell contains b'*alpha which is useful to compute
%               duality gap, etc.
%   progress    Primal objective, duality gap etc as the algorithm progresses,
%               can be used to visualize the convergence.
%
% Authors: Simon Lacoste-Julien, Martin Jaggi, Mark Schmidt, Patrick Pletscher
% Web: https://github.com/ppletscher/BCFWstruct
%
% Relevant Publication:
%       S. Lacoste-Julien, M. Jaggi, M. Schmidt, P. Pletscher,
%       Block-Coordinate Frank-Wolfe Optimization for Structural SVMs,
%       International Conference on Machine Learning, 2013.

% == getting the problem description:
phi = param.featureFn; % for \phi(x,y) feature mapping
loss = param.lossFn; % for L(ytruth, ypred) loss function
if isfield(param, 'constraintFn')
    % for backward compatibility with svm-struct-learn
    maxOracle = param.constraintFn;
else
    maxOracle = param.oracleFn; % loss-augmented decoding function
end

patterns = param.patterns; % {x_i} cell array
labels = param.labels; % {y_i} cell array
n = length(patterns); % number of training examples

% == parse the options
options_default = defaultOptions(n);
if (nargin >= 2)
    options = processOptions(options, options_default);
else
    options = options_default;
end

% general initializations
lambda = options.lambda;
beta = options.beta;
phi1 = phi(param, patterns{1}, labels{1}); % use first example to determine dimension
d = length(phi1); % dimension of feature mapping
using_sparse_features = issparse(phi1);
progress = [];

% === Initialization ===
% set w to zero vector
% (corresponds to setting all the mass of each dual variable block \alpha_(i)
% on the true label y_i coordinate [i.e. \alpha_i(y) =
% Kronecker-delta(y,y_i)] using notation from Appendix E of paper).

% w: d x 1: store the current parameter iterate
% wMat: d x n matrix, wMat(:,i) is storing w_i (in Alg. 4 notation) for example i.
%    Using implicit dual variable notation, we would have w_i = A \alpha_[i]
%    -- see section 5, "application to the Structural SVM"
if using_sparse_features
    model.w = sparse(d,1);
    model.w_y = sparse(d,1);
    wMat = sparse(d,n); 
else
    model.w = zeros(d,1);
    model.w_y = zeros(d,1);
    wMat = zeros(d,n); 
end

model.ell = 0; % this is \ell in the paper. Here it is assumed that at the true label, the loss is zero
         % Implicitly, we have ell = b' \alpha

% logging
progress.primal = [];
progress.dual = [];
progress.gap = [];
progress.eff_pass = [];
progress.train_error = [];
if (isstruct(options.test_data) && isfield(options.test_data, 'patterns'))
    progress.test_error = [];
end

fprintf('running batch FW on %d examples. The options are as follows:\n', length(patterns));
options

tic();

%list of active vertices
model.w_v = [];
model.w_y_v= [];
model.ell_v = [];
% === Main loop ====
for p=1:options.num_passes
    k = p-1; % same k as in paper
    
    if using_sparse_features
        w_s = sparse(d,1);
    else
        w_s = zeros(d,1);
    end
    ell_s = 0;
    %Constructing FW step
    for i = 1:n
        % solve the loss-augmented inference for point i
        ystar_i = maxOracle(param, model, patterns{i}, labels{i});        
        % define the update quantities:
        % [note that lambda*w_s is subgradient of 1/n*H_i(w) ]
        % psi_i(y) := phi(x_i,y_i) - phi(x_i, y)
        psi_i =   phi(param, patterns{i}, labels{i}) ...
                - phi(param, patterns{i}, ystar_i);
        w_s = w_s + 1/n * psi_i;
        loss_i = loss(param, labels{i}, ystar_i);
        ell_s = ell_s + 1/n*loss_i;
        % sanity check, if this assertion fails, probably there is a bug in the
        % maxOracle or in the featuremap
        assert((loss_i - model.w'*psi_i) >= -1e-12);
    end
    %w_y = Ay_s
    w_y_s = w_s;
    % t he true w_s only contains the opposite of the biggest (in absolute value) coordinate times beta 
    w_y = model.w_y;
    [~,i] = max(abs(w_y));
    w_s(1:end ~= i) = 0;
    w_s(i) = beta*sign(w_y(i));
    assert(sum(abs(w_s))<=beta)

    %Constructing Away step
    if ~isempty(S_t) 
        id_A   = away_step(Fz, S_t, I_active);
        v_A    = S_t(:, id_A);
        d_A    = z - v_A;
        alpha_max = alpha_t(id_A);
    else
        fprintf('error: empty support set at step (it=%f)\n', it);
    end


    % construct direction (between towards and away):
    
    if isempty(S_t) || - gap <= d_A' * grad
        is_aw = false; % this is a pure FW step
        d = d_FW; 
        max_step = 1;
    else
        is_aw = true; % this is an away step
        number_away = number_away+1;
        d = d_A;
        max_step = alpha_max / (1 - alpha_max);
    end
    
    gap_AFW = -Fz'*d;
    if safety
        gamma_max = 2/(k+2);
    else
        gamma_max = 1;
    end
    gamma = gamma_coef*gap_AFW;
    gamma_k = min([gamma,gamma_max,max_step]);
    
    step = gamma_k;

    % doing steps and updating active set:
    
    if is_aw
        % away step:
        alpha_t = (1+step)*alpha_t; % note that inactive should stay at 0;
        if abs(step - max_step) < 10*eps
            % drop step:
            number_drop = number_drop+1;
            alpha_t(id_A) = 0;
            I_active(I_active == id_A) = []; % remove from active set
            %TODO: could possibly also remove it from S_t
        else
            alpha_t(id_A) = alpha_t(id_A) - step;
        end
    else
        % FW step:
        alpha_t = (1-step)*alpha_t;
        
        % is this s_FW a new vertex?
        h = hashing(s_FW);
        if ~mapping.isKey(h)
            % we need to add vertex in S_t:
            max_index = max_index + 1;
            mapping(h) = max_index;
            S_t(:,max_index) = s_FW;
            id_FW = max_index;
            alpha_t(id_FW) = step; % this increase size of alpha_t btw
            I_active = [I_active, id_FW];
        else
            id_FW = mapping(h);
            if alpha_t(id_FW) < eps
                % we already had atom in 'correction poytope', but it was not
                % active, so now track it as active:
                I_active = [I_active, id_FW];
            end
            alpha_t(id_FW) = alpha_t(id_FW) + step;
        end
        
        % exceptional case: stepsize of 1, this collapses the active set!
        if step > 1-eps;
            I_active = [id_FW];
        end
    end

    % compute duality gap:
    
    gap_omega = (model.w-w_s)'*(-model.w_y);
    gap_y = ell_s -  model.ell -  model.w'* (w_y_s - model.w_y);
    gap = gap_omega + gap_y;    % keyboard
    % get the step-size gamma:
    %%%% !!!! No line search for saddle points
    if (options.do_line_search) 
        % analytic line-search for the best stepsize [by default]
        % formula from Alg. 2:
        gamma_opt = gap / (lambda*( (model.w - w_s)'*(model.w - w_s) +eps));
        % +eps is to avoid division by zero...
        gamma = max(0,min(1,gamma_opt)); % truncate on [0,1]
    else
        % we use the fixed step-size schedule
        gamma = 2/(k+2);
    end
    
    % stop if duality gap is below threshold:
    if gap <= options.gap_threshold
        fprintf('Duality gap below threshold -- stopping!\n')
        fprintf('current gap: %g, gap_threshold: %g\n', gap, options.gap_threshold)
        fprintf('Reached at iteration %d.\n', k)
        break % exit loop!
    else
        fprintf('Duality gap check: gap = %g at iteration %d\n', gap, k)
    end
    
    % finally update the weights and ell variables
    model.w_y = (1-gamma)*model.w_y   + gamma*w_y_s;
    model.w   = (1-gamma)*model.w   + gamma*w_s;
    model.ell = (1-gamma)*model.ell + gamma*ell_s;
        
    % debug: compute objective and duality gap. do not use this flag for
    % timing the optimization, since it is very costly!
    if (options.debug)
        progress.gap = [progress.gap; gap];
        % primal = -objective_function(model.w, model.w_y,model.ell);
        % % gap = duality_gap(param, maxOracle, model, lambda);
        % % primal = f+gap; % a cheaper alternative to get the primal value
        % train_error = average_loss(param, maxOracle, model);
        % fprintf('pass %d (iteration %d), SVM primal = %f, duality gap = %f, train_error = %f \n', ...
        %     p, k+1, primal, gap, train_error);
        
        % progress.primal = [progress.primal; primal];
        % % progress.dual = [progress.dual; f];
        progress.eff_pass = [progress.eff_pass; k];
        % progress.train_error = [progress.train_error; train_error];
        % if (isstruct(options.test_data) && isfield(options.test_data, 'patterns'))
        %     param_debug = param;
        %     param_debug.patterns = options.test_data.patterns;
        %     param_debug.labels = options.test_data.labels;
        %     test_error = average_loss(param_debug, maxOracle, model);
        %     progress.test_error = [progress.test_error; test_error];
        % end
    end
    
    % time-budget exceeded?
    t_elapsed = toc();
    if (t_elapsed/60 > options.time_budget)
        fprintf('time budget exceeded.\n');
        return
    end
end

end % solverFW


function options = defaultOptions(n)

options = [];
options.num_passes = 200;
options.do_line_search = 1;
options.time_budget = inf;
options.debug = 0;
options.lambda = 1/n;
options.test_data = [];
options.gap_threshold = 0.1;

end % defaultOptions
