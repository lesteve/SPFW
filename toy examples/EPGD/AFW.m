function p = AFW(z,z_0,niter,epsilon)
% Project on the unit cube with AFW
%===== ALG CHOICE ========
% use niter == 0 for projection with epsilon precision.
% use niter == k for k step of AFW 


% ===== FUNCTIONS =====
F = @(p) 2*(p-z) ; 
S = @(F) -(sign(F+eps)-1)*.5;

% ==== INITIALIZATION ====
x_0 = z_0;
S_0 = z_0;
alpha_0 = 1;

iter = 1;
tol = epsilon;

p           = x_0;
S_t         = S_0;
alpha_t     = alpha_0;

% each column of S_t is a potential vertex
% I_active -> contains index of active vertices (same as alpha_t > 0)
mapping = containers.Map();
% this will map vertex hash to id in S_t (to see whether already seen vertex)
% alpha_t(i) == 0 implies vertex is not active anymore...
% alpha_t will be a weight vector so that x_t = S_t * alpha_t

% constructing mapping:
max_index = size(S_t,2); % keep track of size of S_t
for index = 1:max_index
    mapping(hashing(S_t(:,index))) = index; 
end
I_active = find(alpha_t > 0);
gap_best = 10^10;

% tracking results:
number_drop = 0; % counting drop steps (max stepsize for away step)
number_away = 0;
if niter == 0
	Kmax = 1/epsilon;
else
	Kmax = niter;
end
% ==== MAIN ALGORITHM ====
for k = 0:Kmax
    Fz = F(p);
    Sz = S(Fz); % FW-corner
    s_FW = Sz;
    grad = Fz; 
    d_FW = s_FW - p;
    gap= - d_FW' * grad;
    if gap < tol
        % fprintf('Converged at iter %d!\n', iter)
        if gap < 0 
            assert false
        end
        break
    end
    % away direction search:
    if ~isempty(S_t) 
        id_A   = away_step(Fz, S_t, I_active);
        v_A    = S_t(:, id_A);
        d_A    = p - v_A;
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
    gamma_line_search = (z - p)'* d / norm(d)^2;
    step = max(min(max_step,gamma_line_search),0);
    % doing steps and updating active set:
    
    if is_aw
        % away step:
        alpha_t = (1+step)*alpha_t; % note that inactive should stay at 0;
        if abs(step - max_step) < 10*eps
            % drop step:
            number_drop = number_drop+1;
            alpha_t(id_A) = 0;
            I_active(I_active == id_A) = []; % remove from active set
        else
            alpha_t(id_A) = alpha_t(id_A) - step;
        end
    else
        % FW step:
        alpha_t = (1-step)*alpha_t;
        
        % is this s_FW a new vertex?
        h = hashing(2*s_FW-1);
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
    p = p + step * d;
    iter = iter+1;
end

