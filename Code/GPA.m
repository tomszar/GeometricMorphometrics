function [shapematrix, cs] = GPA(landmarkmatrix, symmetrize)
  %Returns shape coordinates after applying a Generalized Procrustes Analysis
  
  %GPA runs a Generalized Procrustes Analysis on a set of 3D landmarks
  %The landmark matrix should be a n by k*3 numeric matrix where n is the 
  %sample size and k is the number of landmarks arranged in the form of 
  %x1, y1, z1, x2, y2, z2, ..... , xk, yk , zk
  %Set symmetrize 1 to return a symmetrized shape coordinate data
  
  %The algorithm runs as follows:
  %1. CENTERING AND SCALING: 
  %First, center each configuration at the origin by substracting the 
  %coordinates of the centroid from the corresponding coordinates at each 
  %landmark.
  %Then, scale configuration to unit centroid size,by dividing each coordinate 
  %of each landmark by the centroid size
  
  %2. ROTATION
  %Rotates each configuration w.r.t the first one using a SVD approach.
  %A consensus is computed and all configurations are rotated with respect to it
  %New iterations of the procedure are applied until the difference in the 
  %residual sum of squares is less than 10^-4
  
  form = reshape_coords(landmarkmatrix);
  
  if nargin == 1
    symmetrize = 0;
  endif
  
  %Create form matrices to ease computation
  p    = size(landmarkmatrix,2)/3; %points
  k    = 3; %dimensions 
  n    = size(landmarkmatrix,1); %sample
  
  %GPA proper procedure

  %%%% 1. CENTERING AND SCALING%%%%
 
  %First, center each configuration at the origin
  %by substracting the coordinates of the centroid from the corresponding 
  %coordinates at each landmark_ids
  %Then, scale configuration to unit centroid size,
  %by dividing each coordinate of each landmark by the centroid size

  %Getting the centroid and centering matrix
  %form_i = zeros(size(form));
    
  if symmetrize == 1
    formr = form;
    formr(:,1,:) = formr(:,1,:) * -1;
    %flipping the x coordinates to maintain correspondence of landmarks
    %formr(:,1,:) = flipud(formr(:,1,:));
    form  = cat(3, form, formr); 
    n     = size(form, 3);
  endif
  
  cs     = zeros([1 n]);
  
  '1. Centering and scaling....'
  centroid = [mean(form(:,1,:)), mean(form(:,2,:)), mean(form(:,3,:))];
  form     = form - centroid;
  cs(1,:)  = sqrt(sum(reshape(form, [1 p*k n]).^2));
  form     = reshape_coords(form);
  form     = form ./ cs';
  form     = reshape_coords(form);
   
  %%%% 2.ROTATION %%%%
  '2. Rotating...'   
  Y = form(:,:,1);
  X = form(:,:,2:n);
  
  '    - W.r.t. the first configuration...'
  Z = cellfun(@(x) Y'*x, num2cell(X,[1 2]),'UniformOutput',false);
  Z = cat(3,Z{:});
  [U,~,V] = cellfun(@(x) svd(x), num2cell(Z,[1 2]),'UniformOutput',false);
  for i = 2:n
    form(:,:,i) = X(:,:,i-1) * (cell2mat(V(:,:,i-1)) * cell2mat(U(:,:,i-1))');
  endfor
    
  '    - W.r.t. the consensus...'
  %Rotations with respect to the mean
  fit = 1; %Setting initial fit to enter the while loop
  ps  = ones([1 n]);
  ps_dot = ones([1 n]);
  Y  = sum(form, 3) / n; %Setting the new reference to the average shape
  sr = n*(1 - trace(Y*Y') ); %Initial residual SS

  while fit > 10^-4
    '    - Iterate through...'
    X  = form(:,:,:);
    Z = cellfun(@(x) Y'*x, num2cell(X,[1 2]),'UniformOutput',false);
    Z = cat(3,Z{:});
    [U,~,V] = cellfun(@(x) svd(x), num2cell(Z,[1 2]),'UniformOutput',false);
    
    for i = 1:n
      form(:,:,i) = ps(1,i) * X(:,:,i) * (cell2mat(V(:,:,i)) * cell2mat(U(:,:,i))');
    endfor
  
    %New consensus
    %Faster sepparating the computations in a and b
    %%doing trace(X * Y') with sum(sum(X.*conj(Y),2))
    Ydot   = sum(form, 3) / n;
    c      = sum(sum(Ydot.*conj(Ydot),2));
    a      = cellfun(@(x) sum(sum(x.*conj(Ydot),2)), num2cell(form,[1 2]),'UniformOutput',false);
    b      = cellfun(@(x) sum(sum(x.*conj(x),2)) * c, num2cell(form,[1 2]),'UniformOutput',false);
    lefts  = sqrt(cell2mat(a) ./ cell2mat(b));    
    ps_dot = reshape(lefts, [size(ps)]) .* ps;
    for i = 1:n
      form(:,:,i) = lefts(i) * form(:,:,i);
    endfor
  
    %New consensus and comparison between residuals
    Ydot2 = sum(form, 3) / n;
    sr2 = n*(1 - trace(Ydot2*Ydot2') );
    Y = Ydot2;
    ps  = ps_dot;
    fit = abs(sr - sr2);
    sr  = sr2;
    fit
    
  endwhile 

  if symmetrize == 1
    shapematrix = zeros(size(form(:,:,1:(n/2) ) ) );
    %Removing repeated cs estimations
    cs = cs(1:(n/2));
    for i = 1:(n/2)
      shapematrix(:,:,i) = ( form(:,:,i) + form(:,:, i + (n/2)) ) /2;
    endfor 
  elseif symmetrize == 0
    shapematrix = form;
  endif
  
  shapematrix = reshape_coords(shapematrix);
   
endfunction

function [ new_coords ] = reshape_coords( coords )
  %Reshape 3D landmark coordinate data
  
  %Reshape_coords takes a p x 3 x n array, or a n x p matrix and returns
  %the other one, where p is the number of landmarks, and n is the sample size
  
  k = 3;
  
  if size(coords, 3) > 1
    p = size(coords, 1);
    n = size(coords, 3);
    new_coords = zeros(n, p*3);
    
    for i = 1:n
      matrix          = coords(:,:,i);
      new_coords(i,:) = reshape(matrix.', 1, []);
    endfor
    
  elseif size(coords, 3) == 1
    p = size(coords, 2)/3;
    n = size(coords, 1);
    new_coords = zeros(p, k , n);
    
    for i = 1:n
      new_coords(:,:,i) = reshape(coords(i,:), [k, p] )';
    endfor
    
  endif
  
endfunction
