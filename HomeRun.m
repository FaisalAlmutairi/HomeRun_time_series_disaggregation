function [s, cost] = HomeRun(A, H, y, rho, lambda)
% This code solve HomeRun formulation via ADMM.
% Details are in the paper:
% "HomeRun: Scalable Sparse-Spectrum Reconstruction of Aggregated Historical Data"
%
% This function solves the following problem via ADMM:
%
%          minimize     ||s||_1 + (lbd/2)||H*iDCT(s)||_2
%          subject to   A*iDCT(s) = y
%                       iDCT(s) > = 0
%
% iDCT is inverse Discrete Cosine Transform (DCT). DCT and iDCT are performed 
% using FFT using mirt_dctn and mirt_idctn functions written by Andriy Myronenko.         
%
% The solution is returned in s. rho is the augmented Lagrangian parameter.
%
% Faisal Almutairi, almut012@umn.edu , Jan. 2018.


A = sparse(A);
H = sparse(H);
y = sparse(y);

% some constants
MAX_ITER = 3000;
cost = zeros(1,MAX_ITER);
[m, n] = size(A);

% initial values 
s = sparse(n,1);
r = sparse(n,1);
u1 = sparse(m,1);
u2 = sparse(n,1);
u3 = sparse(n,1);

% compute Cholesky decomposit
Z = (rho*(A'*A) + 2*rho*speye(n) + lambda*(H'*H));
L = chol(Z,'lower');

for iter = 1:MAX_ITER 
    
    %%ADMM Steps (z and r are  auxiliary variables)
    %% Step 1: update z (least square closed form)
    bs = rho*A'*(y - u1) + rho*(r - u2) + rho*sparse(mirt_idctn(full(s - u3)));
    t = L'\(L\bs);
    z = sparse(mirt_dctn(full(t)));
    
    %% Step 2: update r and s      
    % r-update (non-negativity projection)
    r = t + u2;
    r(r<0) = 0;
    % s-update (soft thresholding)
    sold = s;
    s = sparse(max(0,(z + u3)-(1/rho)) - max(0, -(z + u3)-(1/rho)));        
    %% Step 3: update the dual variables
    u1 = u1 + A*t - y;
    u2 = u2 + t - r;
    u3 = u3 + (z - s);
    
    %% calculate cost 
    cost(iter) = objective(s,H,lambda);
    costold = objective(sold,H,lambda);
    if abs(cost(iter)-costold)/cost(iter) < eps
        break
    end
 
end

end


function obj = objective(s,Hs,lbd)
    obj = norm(s,1) +(lbd/2)*(norm(Hs*mirt_idctn(full(s)),2))^2;
end