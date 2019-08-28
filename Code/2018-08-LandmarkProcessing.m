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

%Remove the string 'PSU' from the landmark_ids file, and remove the starting 0
landmark_ids = erase(landmark_ids, 'PSU');
landmark_ids = regexprep(landmark_ids,'^0*','');

%Exporting total set of IDs
cd(folder.save);
cell2csv("total_faceshape_ids.txt", landmark_ids)

%Loading SNP IDs after combination of datasets
cd(folder.databases);
geno_ids = textread ('SNP_2729_ids.txt', '%s', 'delimiter' , ' ');

%Retain IDs from genotype data into shape analysis
%Then, set the intersection
[~ , keep] = intersect(landmark_ids, geno_ids);

%Keep those only in the landmark_matrix
landmark_matrix = landmark_matrix(keep,:);
landmark_ids    = landmark_ids(keep,:);

%Run GPA
[shape_matrix, cs] = GPA(landmark_matrix);

%Read Covariates and generate the average male and female face from the shape_matrix
[ID, sex, ~, ~, ~] = textread ('Covariates.csv', '%s %s %f %f %f', 'delimiter' , ',', 'headerlines', 1);
[~ , keep] = intersect(ID,landmark_ids);

sex = sex(keep,:);
ID  = ID(keep,:);

female_avg = mean(shape_matrix(strcmp(sex, 'Female'),:) );
male_avg   = mean(shape_matrix(strcmp(sex, 'Male'),:) );
total_sex_avg_shape = [female_avg; male_avg];

%Delete landmark_matrix
clear landmark_matrix;

%Run fastPCA on whole sample
[V,S,score,eigenvals,percent] = fastPCA(shape_matrix, 100);

%Retain PCs with eigenvalues higher than 1%
index = find(percent > 1);

%plot(eigenvals,'bo-');
%hold on;
%plot(mean_eigenvals,'ro-');
%hold off;

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
csvwrite("total_avg_sex_shape.txt", total_sex_avg_shape)
csvwrite("means.txt", mu')
csvwrite("cs.txt", cs') 
csvwrite("scores.csv", score) 
csvwrite("eigenvalues.csv", [eigenvals, percent'] ) 
csvwrite("eigenvectors.csv", V )
csvwrite("facets.csv", landmark_facets ) 
