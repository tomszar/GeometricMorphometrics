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
%Loading common IDs from combination of datasets
geno_ids = textread ('common_ids.txt', '%s', 'delimiter' , ' ');

%Retain IDs from genotype data into shape analysis
%First, remove the string 'PSU' from the landmark_ids file, and remove the starting 0
landmark_ids = erase(landmark_ids, 'PSU');
landmark_ids = regexprep(landmark_ids,'^0*','');
%Then, set the intersection
[~ , keep] = intersect(landmark_ids, geno_ids);

%Keep those only in the landmark_matrix
landmark_matrix = landmark_matrix(keep,:);
landmark_ids    = landmark_ids(keep,:);

%Run GPA
[shape_matrix, cs] = GPA(landmark_matrix);

%Delete landmark_matrix
clear landmark_matrix;

%Run fastPCA on whole sample, with PA and 999 runs for PA
[V,S,score,eigenvals,percent,cutoff_eigenvals,index] = fastPCA(shape_matrix, 100, 'PA', 999);

%Looking at the number of relevant PCs
plot(eigenvals,'bo-');
hold on;
plot(mean_eigenvals,'ro-');
hold off;

%Retain only relevant PCs from PA
score     = score(:,index);
V         = V(:,index);
%eigenvals = eigenvals(index,:);
%percent   = percent(:,index);

%Get the mean to export
mu = mean(shape_matrix);

%Saving files
cd(folder.save)
cell2csv("landmark_ids.txt", landmark_ids)
csvwrite("means.txt", mu')
csvwrite("cs.txt", cs') 
csvwrite("scores.csv", score) 
csvwrite("eigenvalues.csv", [eigenvals, percent'] ) 
csvwrite("eigenvectors.csv", V )
csvwrite("facets.csv", landmark_facets ) 
