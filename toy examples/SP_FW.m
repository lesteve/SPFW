% Frank-Wolfe algorithm **with away steps**
%===== ALG CHOICE ========
away = 1
    % use 0 for regular FW
    % use 1 for away step
% saddle point function:
% L(x,y) = x^2 - 2*alpha*x*y - y^2
% ==== STEP SIZE CHOICE =====
% adaptive = 1; % flag to choose adaptive step-size
    % use adaptive == 1 for gap/C
    % use adpative == 2 for line-search...
    % use adaptive == 3 for 2/(2+k)
safety = 0; % flag to use the 2/(k+2) safety when using big step-size...
% ==========================

x = -1:0.05:1;
[X,Y] = meshgrid(x,x);
%L = inline('X.^2 - 2*c*X.*Y - Y.^2', 'X','Y','c');

% dim= 30;
rng(2)
% cst = .1
M = cst.* (2.*rand(dim,dim)-1); % for 

% Frank-Wolfe; variational inequality...

% original function:
%F = @(x,y,c) [2*x - 2*c*y; 2*c*x+2*y]; % F(z)
% above converges at 1/k easily...

% Now considering translated version:
%L = inline('0.5*(X-a).^2 - c*(X-a).*(Y-b) - 0.5*(Y-b).^2', 'X','Y','c');
% put global saddle point at (a,b), not sure now what will happen!
metaL = @(x,y,M,a,b,strong) 0.5.*strong.*sum((x-a).^2) - (x-a)'*M*(y-b)-0.5.*strong.*sum((y-b).^2); 
metaF = @(x,y,M,a,b,strong) [strong.*(x-a) - M*(y-b); M'*(x-a)+strong.*(y-b)]; 
% a = 2*ones(dim,1)-1; b=2*ones(dim,1)-1; 
% strong=8.4; % strong = convexity for each term
%100 50 35 31.8
F = @(x,y,M) metaF(x,y,M,a,b,strong); % F(z)
L = @(x,y,M) metaL(x,y,M,a,b,strong);

%Z = L(X,Y,c);
%surf(X,Y,Z)

%
% [m,i] = max(abs(F)
S = @(F) -sign(F+eps);

% gap = <S(z) - z, -F(z) >


%C = strong*4; %L*diam^2 -- old constant
mu = strong; % strong convexity constant

% distance_x = min(abs(a-1), abs(a+1));
% distance_y = min(abs(b-1), abs(b+1));
% distance_boundary = min(distance_x,distance_y);
% mu_FW = mu*(distance_boundary^2);
% M_L = 2*c*2; % L_xy D_x D_y
% alpha = sqrt(2)*M_L / mu_FW; % alpha in NIPS Gauthier submission

% ** new correct constants **
% C = C_x + C_y + M_xy * M_yx [D_x^2/mu_x + D_y^2/mu_y]
% note diameter here is 2 for each dimension...

%
delta = min(norm(min([a - 1, a + 1]'),2),norm(min([b-1,b+1]),2));
Pwidth = 1/sqrt(dim);
diam = 2* sqrt(dim);
if away==1
    alpha = 1/2 - 1/(strong*Pwidth^2)*max(max(M))*diam/(strong)
else 
    alpha = 1 - 1/(strong*delta^2)*max(max(M))*diam/(strong)
end
% assert(alpha >0)
% Maybe the "good" C
if adaptive ==4
    C = 2*strong*diam^2+ max(max(abs(M)))^2*(diam^2/strong+diam^2/strong)
else 
    C = strong*diam^2/alpha;
end

if adaptive > 0
    gamma_coef = 1/(C)   
end





% x_0 = [zeros(dim,1); zeros(dim,1)];
% x_0(1)=1;
% x_0(dim+1)=1;
x_0 = [sign(2*rand(dim,1)-1); sign(2*rand(dim,1)-1)];
S_0 = x_0;
alpha_0 = 1;

iter = 1;
tol = 1e-5;
G = [];

z           = x_0;
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

for k = 0:Kmax
    
    Fz = F(z(1:dim),z(dim+1:end),M);
    

    Sz = S(Fz); % FW-corner
    % [m,i_x] = max(abs(Fz(1:dim)));
    % [m,i_y] = max(abs(Fz(dim+1:end)));
    % Sz = zeros(2*dim,1);
    % Sz(i_x) = -sign(Fz(i_x));
    % Sz(i_y+dim) = -sign(Fz(i_y+dim));

    s_FW = Sz;
    grad = Fz;
    
    d_FW     = s_FW - z;
    gap= - d_FW' * grad;
    gap_best = min(gap,gap_best);
    G(iter) = gap_best; % keep track of FW gaps..
    
    if gap < tol
        fprintf('Converged at iter %d!\n', iter)
        break
    end
    if away 
         % away direction search:

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
        if adaptive ==3
            gamma_k = 2./(2+iter);
        end
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
    else
            step = min(1,gamma_coef*gap);
            d = d_FW;
    end

    z = z + step * d;
    % assert(max(abs(z(1:dim)))<=1 +1000*eps)
    % assert(max(abs(z(dim:end)))<=1 +1000*eps)
    %z = (1-gamma_k)*z + gamma_k * Sz;

    iter = iter+1;
end

% rng(C/100)
% if adaptive
%     semilogy(G,'color',[rand,rand,rand])
% else
%     semilogy(G,'color',[rand,rand,rand])
%     %loglog(G)
% end



    


