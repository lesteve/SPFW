function [G,iter] = EPGD(Kmax,mu,M,cst,dim,a,b,alg,adaptive)
% Projected Extra-gradient Descent
% Solve the saddle point problem 
% L(x,y) = mu*\|x\|^2 - x^\top M y - mu*\|y\|^2
% Using the Projected Extra-gradient Descent see Korpelevich 197
%
% Each step does a gradient descent step (for SP).
% Inputs:
%   Kmax: The number of iteration
%   mu: The value of the mu convexity constant     
%   cst: The matrix M \sim Unif([-cst,cst]^2*dim)
%   dim: Dimension of the problem
%   a: position of the saddle point in \X
%   b: position of the saddle point in \Y
%   away: chose if the algorithm is SP-FW of SP-AFW
%   adaptive: define the step size of the algorithm
%       
% %===== ALG CHOICE ========
    % use alg == 0 for exact projection
    % use alg == 1 for 1 step of AFW for inexact projection
    % use alg == 2 for epsilon/2 precision of inexact projection with AFW
% ==== STEP SIZE CHOICE =====
    % use adaptive == 0 for 1/L
    % use adaptive == 1 for 1/sqrt(1+k)
    % use adaptive == 2 adaptive step size
if alg == 0
    proj = @(z,z_0) exact_proj(z);
elseif alg == 1
    proj = @(z,z_0) AFW(z,z_0,1,1e-10);
else
    proj = @(z,z_0) AFW(z,z_0,0,epsilon);
end  

if adaptive == 0
    step = @(k) .1/mu;
elseif adaptive == 1
    step = @(k) 1/(max(norm(M),mu)*sqrt(k+1));
else
    step = @(k) k; %see later
end  




% Definition of the LMO 
metaF = @(x,y,M,a,b,mu) [mu.*(x-a) - M*(y-b); M'*(x-a)+mu.*(y-b)]; 
F = @(x,y) metaF(x,y,M,a,b,mu); % F(z)
S = @(F) -(sign(F+eps)-1)*.5;

% ==== INITIALIZATION ====
x_0 = [rand(2*dim,1)];
S_0 = x_0;
alpha_0 = 1;

iter = 1;
tol = 1e-10;
G = [];

z           = x_0;
S_t         = S_0;
alpha_t     = alpha_0;


% ==== MAIN ALGORITHM ====
for k = 0:Kmax
    Fz = F(z(1:dim),z(dim+1:end));
    Sz = S(Fz); % FW-corner to compute the gap
    grad = Fz; 
    d_FW     = Sz - z;
    gap= - d_FW' * grad ;
    if iter == 1 
        gap_best = gap;
    end
    if gap < tol
        fprintf('Converged at iter %d!\n', iter)
        if gap < 0 
            assert false
        end
        break
    end
    gap_best = min(gap_best,gap);
    G(iter) = gap_best; % keep track of FW gaps..
    z_bar = proj(z - step(iter) * Fz,z);
    Fz_bar = F(z_bar(1:dim),z_bar(dim+1:end));
    z = proj(z - step(iter) * Fz_bar,z);
    iter = iter +1;
    % keyboard
end


    


