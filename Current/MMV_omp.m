function [S] = MMV_omp(A,B,K)

% OMP alogrithm for MMV: B=AX
% MMV-OMP:Orthogonal Matching Pursuit Algorithm Implementation
%
% Inputs
%   A           : Sensing Matrix (size MxN)
%   B           : Measurement Vector (size MxL) B = AX
%   K           : sparsity level
% Outputs
%   Supp        : the support of X
% Purpose       : Implement Orhtogonal Matching Pursuit Algorithm 

[mA,nA] = size(A);
[mY,nY] = size(B);

residual = B;
Supp = [];
iter = 1;

NormACols = sqrt(diag(A'*A));

while ((iter <= K))
    % Matching step
    Z_1 = A'*residual;
    Z = sqrt(sum(abs(Z_1).^2,2))./NormACols;

    [maxVal, maxPos] = max(Z);
    BestLoc = maxPos(1);
    
    % Update support
    Supp = [Supp BestLoc];
    As = A(:,Supp);

    solution = As*(As\B);
    residual = B-solution;

    % increment iter number
    iter = iter+1;
end

% Construct solution
% Xr = zeros(nA,nY);
% Xr(Supp,:) = pinv(As)*B;

S=pinv(As)*B;
