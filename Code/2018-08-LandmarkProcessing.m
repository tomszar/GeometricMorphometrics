% This script will run a Generalized Procrustes Analysis from the GPA function
% and a PCA from the fastPCA function
% First, we will keep only those individuals with genotype data, from the adapt_geno.fam file

%Adding code from current directory
addpath (pwd())

%Load pkg
pkg load statistics
pkg load io
pkg load geometry

%Get folders and load databases
cd .. 
folder.root      = pwd();
folder.databases = strcat(folder.root, '/DataBases');
folder.save      = strcat(folder.root, '/Results');
cd(folder.databases);
load('landmark_matrix.mat');
load('landmark_ids.mat');
load('landmark_facets.mat');
%Loading IDs
[~, geno_ids, ~,~,~,~] = textread ('adapt_geno.fam', '%s %s %f %f %f %f', 'delimiter' , ' ');
%ids    = ids(2:end,:);
%groups = groups(2:end,:);

%Retain IDs from genotype data into shape analysis
%First, remove the string 'PSU' from the landmark_ids file, and remove the starting 0
landmark_ids = erase(landmark_ids, 'PSU');
landmark_ids = regexprep(landmark_ids,'^0*','');
%Then, set the intersection
[~ , keep] = intersect(landmark_ids, geno_ids);

%Keep those only in the landmark_matrix
landmark_matrix = landmark_matrix(keep,:);
landmark_ids    = landmark_ids(keep,:);


%PCA (Run this with random equal samples from each group)
%Count group members
%[C,ia,ic] = unique(groups);
%group_counts  = accumarray(ic,1);
%group_counts

%Get random sample of 120 per group and run PCA on that
%row = 1
%sample = 120
%random_sample = zeros(size(C)(1) * sample, size(shape_matrix )(2));
%for i = 0:5
%  reduced_sample = shape_matrix(groups == i,:);
%  start = 1 + (i *sample);
%  fin   = row * sample;
%  random_sample(start:fin, :) = reduced_sample(randperm(sample), :);
%  row = row + 1;
%endfor
%Run PCA on random sample
[V,S,score,mean_eigenvals,index] = fastPCA(shape_matrix, 100, 1);
%Get score to whole dataset
mu = mean(shape_matrix);
Xm =  bsxfun(@minus, shape_matrix, mu);
V     = V(:,1:index);
total_scores = Xm * V;

%Run PCA with whole sample
%[V,S,score,mean_eigenvals,index] = fastPCA(shape_matrix, 100, 1);
#Looking at the number of relevant PCs
plot(diag(S),'bo-');
hold on;
plot(mean_eigenvals,'ro-');
hold off;
%Getting PCA scores for relevant PCs
%mu    = mean(shape_matrix);
%score = score(:,1:index);
%V     = V(:,1:index);
%Eigenvalue is (diag(S).^2)/(n-1)

%Testing accuracy of PCA
%Calculating euclidean distance between original and PCA transformed faces
fromPCA   = total_scores * V' + mu;
distances = distancePoints(shape_matrix, fromPCA, 'diag');
distances_rand = distancePoints(shape_matrix, fromPCA(randperm(size(fromPCA,1)), :), 'diag');
hist(distances,'bo-');
hold on;
hist(distances_rand,'ro-');
hold off;

#Plot average faces
avg_scores = zeros(7, size(score,2));
for (i = 1:7)
  [~ , keep] = intersect(landmark_ids, ids(groups==i));
  avg_scores(i,:) = mean(score(keep,:));
endfor
from_avgscores = avg_scores * V' + mu;

p = size(fromPCA, 2)/3;
coords = reshape(from_avgscores(5,:), [3, p] )';
position = 0;
dist = [0.8 0.8 0.8];
trisurf(landmark_facets, coords(:,1), coords(:,2), coords(:,3), 'facecolor', dist);
view(2);
axis off;

%Saving files
cd(folder.save)
cell2csv("landmark_ids.txt", landmark_ids)
csvwrite("means.txt", mu')
csvwrite("cs.txt", cs') 
csvwrite("scores.csv", total_scores) 
csvwrite("eigenvalues.csv", diag(S) ) 
csvwrite("eigenvectors.csv", V ) 
csvwrite("facets.csv", landmark_facets ) 
