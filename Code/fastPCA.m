function [V,S,scores,eigenvals,percent,mean_eigenvals,index] = fastPCA(X, K, Panalysis, runs = 99)
% Runs a fastPCA by using a randomized svd algorithm.
% The X matrix is first centered and then a rsvd is applied
% A parallel analysis can be done to determine the number of PCs to retain
% 
% Usage: 
%   Input:  
%   * X         : n by p matrix to obtain the PCA from
%   * K         : number of components to keep
%   * Panalysis : wheter to run a parallel analysis or not (0/1)
%   * runs      : number of replications to obtain the eigenvalues from (default 99)
%
%   Output:
%   * V              : right singular matrix
%   * S              : singular values
%   * scores         : PCA scores from U*S
%   * eigenvals      : eigenvalues calculated from (diag(S).^2)/(n-1)
%   * percent        : percentage of variance explained by each eigenvalue, 
%                      empirically calculated   
%   * mean_eigenvals : if Panalysis == 1, mean_eigenvals is the average eigenvalues 
%                      from a PA analysis of 999 runs
%   * index          : if Panalysis == 1, index is the index of the last eigenvalues 
%                      that is higher than the average one taken from the PA analysis
%

if ~exist('Panalysis', 'var')
  Panalysis = 0;
endif

%Getting n and center matrix
n  = size(X)(1);
mu = mean(X);
Xm = bsxfun(@minus, X, mu);

#Calculate svd from matrix, and generate PCA scores, eigenvals, and % explained
[U, S, V] = rsvd(Xm, K); 
scores    = U*S;
eigenvals = (diag(S).^2)/(n-1);
percent   = var(scores)  / sum(var(Xm)) * 100;

%Wheter to run the parallel analysis or not
if (Panalysis == 1)
  [mean_eigenvals, index] = PA(Xm, K, runs);
elseif (Panalysis == 0)
  mean_eigenvals = NA;
  index = NA;
endif

endfunction


function [U,S,V] = rsvd(A,K)
%-------------------------------------------------------------------------------------
% Random SVD
% Extremely fast computation of the truncated Singular Value Decomposition, using
% randomized algorithms as described in Halko et al. 'finding structure with randomness
%
% usage : 
%
%  input:
%  * A : matrix whose SVD we want
%  * K : number of components to keep
%
%  output:
%  * U,S,V : classical output as the builtin svd matlab function
%-------------------------------------------------------------------------------------
% Antoine Liutkus  (c) Inria 2014
% Modified to Octave by Tomas Gonzalez

[M,N] = size(A);
P = min(2*K,N);
X = randn(N,P);
Y = A*X;
W1 = orth(Y);
B  = W1'*A;
[W2,S,V] = svd(B,'econ');
U = W1*W2;
K = min(K,size(U,2));
U = U(:,1:K);
S = S(1:K,1:K);
V = V(:,1:K);

endfunction

function [mean_eigenvals, index] = PA(X, K, runs)
% Parallel Analysis (PA) to determine the number of components to retain
% The components will be retained if the eigenvalues are greater than their 
% respective eigenvalues from the random data.
%
% Usage:
%   Input:
%   * X    : n by p matrix to return the PA analysis from
%   * K    : number of components to retain
%   * runs : number of replications to obtain the eigenvalues from
%
%   Output:
%   * mean_eigenvals : mean eigenvalues from the random data
%   * index          : last PC to keep based on PA
%

%Calculating real eigenvals
n      = size(X)(1);
[~, S] = rsvd(X, K);
eigenvals_real = (diag(S).^2)/(n-1);
eigenvals_rand = zeros([runs, K]);
shuffledArray  = X;

%Shuffling matrix, and calculating eigenvals
for i = 1:runs
  for dim = 1:size(X,2)
    shuffledArray(:, dim) = X(randperm(size(X,1)), dim);
  endfor
  [~, S] = rsvd(shuffledArray, K); 
  eig    = (diag(S).^2)/(n-1);
  eigenvals_rand(i,:) = eig;
endfor

mean_eigenvals = mean(eigenvals_rand);
index          = find(eigenvals_real' > mean_eigenvals);
  
endfunction
