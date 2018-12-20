function [V,S,scores,eigenvals,percent,cutoff_eigenvals,index,mu,sd] = fastPCA(X, K, stop = 'NA', runs = 99)
%{ 
  Runs a fastPCA by using a randomized svd algorithm.
  The X matrix is first centered and scaled, and then a rsvd is applied
  A stopping rule can be applied to determine the number of PCs to retain
 
  Usage: 
   Input:  
   * X    : n by p matrix to obtain the PCA from
   * K    : number of components to keep
   * stop : type of stopping rule to determine number of PCs to retains 
            (options = 'PA', 'avgPA')
   * runs : number of replications to obtain the eigenvalues for the stopping rule
            (default 99)

   Output:
   * V              : right singular matrix
   * S              : singular values
   * scores         : PCA scores from U*S
   * eigenvals      : eigenvalues calculated from (diag(S).^2)/(n-1)
   * percent        : percentage of variance explained by each eigenvalue, 
                      empirically calculated   
   * cutoff_eigenvals : if stop != 'NA', cutoff_eigenvals is the critical value to
                        determine the retention of PCs
   * index            : if Panalysis == 1, index is the index of the last eigenvalues 
                        that is higher than the average one taken from the PA analysis
   * mu             : list of column means before centering
   * sd             : list of column standard deviations before scaling
%}

%Getting n, and center and scale matrix
n  = size(X)(1);
mu = mean(X);
Xm = bsxfun(@minus, X, mu);
%sd = std(Xm);
%Xm = bsxfun(@rdivide, Xm, sd);

%Getting total variance of matrix
total_var = sum(var(Xm));

%Calculate rsvd from matrix, generate PCA scores, eigenvals, and % explained
[U, S, V] = rsvd(Xm, K); 
scores    = U*S;
eigenvals = (diag(S).^2)/(n-1);
percent   = var(scores)  / total_var * 100;

%Running stopping rule
if strcmp(stop, 'PA')
  [cutoff_eigenvals, index] = PA(Xm, K, runs, type = 'prc');
elseif strcmp(stop, 'avgPA')
  [cutoff_eigenvals, index] = PA(Xm, K, runs, type = 'avg');
elseif strcmp(stop, 'NA')
  cutoff_eigenvals = NA;
  index = NA;
endif

endfunction


function [U,S,V] = rsvd(A,K)
%{
  Random SVD
  Extremely fast computation of the truncated Singular Value Decomposition, using
  randomized algorithms as described in Halko et al. finding structure with randomness
 
  usage : 
   input:
   * A : matrix whose SVD we want
   * K : number of components to keep
 
   output:
   * U,S,V : classical output as the builtin svd matlab function
   
  Antoine Liutkus (c) Inria 2014
  Modified to Octave by Tomas Gonzalez
%}

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

function [cutoff_eigenvals, index] = PA(X, K, runs, type='avg')
%{ 
  Parallel Analysis (PA) to determine the number of components to retain
  The components will be retained if the eigenvalues are greater than their 
  respective eigenvalues from the random data.

  Usage:
    Input:
    * X    : n by p matrix to return the PA analysis from
    * K    : number of components to retain
    * runs : number of replications to obtain the eigenvalues from
    * type : wheter the decision is made from the average eigenval ('avg') or 
             from the 95 percentile ('prc')
 
    Output:
    * cutoff_eigenvals : cutoff eigenvalues from the random data
    * index            : list of PCs to keep based on PA
%}

%Calculating real eigenvals
n      = size(X)(1);
[~, S] = rsvd(X, K);
eigenvals_real = (diag(S).^2)/(n-1);
eigenvals_rand = zeros([runs, K]);

%Getting eigenvals from random data
for i = 1:runs
  rand_matrix = normrnd(0, 1, [size(X)]);
  [~, S] = rsvd(rand_matrix, 100); 
  eigenvals_rand(i,:) = (diag(S).^2)/(n-1);
endfor

%Getting output objects
if (type == 'avg')
  cutoff_eigenvals = mean(eigenvals_rand);  
elseif (type == 'prc')
  cutoff_eigenvals = prctile(eigenvals_rand, 95);
endif

index = find(eigenvals_real' > cutoff_eigenvals);

%{ 
This was for shuffling
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
%}
  
endfunction
