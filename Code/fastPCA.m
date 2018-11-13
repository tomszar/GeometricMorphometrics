function [V,S,scores,mean_eigenvals,index] = fastPCA(X, K, Panalysis)
% Runs a fastPCA by using a randomized svd algorithm.
% The X matrix is first centered and then a rsvd is applied
% 
% Usage: 
%   Input:  
%   X: n by p matrix to obtain the PCA from
%   K: number of components to keep
%   Panalysis: wheter to run a parallel analysis or not (0/1)
%
%   Output:
%   V: right singular matrix
%   S: singular values
%   scores: PCA scores from U*S
%   mean_eigenvals: if Panalysis == 1, mean_eigenvals is the average eigenvalues from a PA analysis of 999 runs
%   index: if Panalysis == 1, index is the index of the last eigenvalues that is higher than the average one taken from the PA analysis
%

if ~exist('Panalysis', 'var')
  Panalysis = 0;
endif

mu = mean(X);
Xm = bsxfun(@minus, X, mu);
%C = cov(Xm);
[U, S, V] = rsvd(Xm, K); 
scores = U*S;

if (Panalysis == 1)
  [mean_eigenvals, index] = PA(Xm, K);
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

function [mean_eigenvals, index] = PA(X, K, runs = 99)
% Parallel Analysis (PA) to determine the number of components to retain
% The components will be retained if the eigenvalues are greater than their respective eigenvalues from the random data.
%
% Usage:
%   Input:
%   X: n by p matrix to return the PA analysis from
%   K: number of components to retain
%   runs: number of replications to obtain the eigenvalues from (default 100)
%
%   Output:
%   mean_eigenvals: mean eigenvalues from the random data
%   index: last PC to keep based on PA
%

[~, S, ~] = rsvd(X, K); 
eigenvals = zeros([runs, K]);
shuffledArray = X;
for i = 1:runs
  for dim = 1:size(X,2)
    shuffledArray(:, dim) = X(randperm(size(X,1)), dim);
  endfor
  [~,S2] = fastPCA(shuffledArray, K);
  eigenvals(i,:) = diag(S2);
endfor
mean_eigenvals = mean(eigenvals);
index = find(diag(S)' > mean_eigenvals)(end);
  
endfunction
