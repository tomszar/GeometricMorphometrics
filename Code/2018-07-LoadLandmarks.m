%Moving to directory
folder.faces = 'C:\Users\tzarzar\Box Sync\PSU_KU_Collaboration\Faces\MAPPED_ALL_1705';
folder.save  = 'C:\Users\tzarzar\Documents';
cd(folder.faces);

%Opening files and saving vertices in new matrix object
files = dir;
dirnames = {files([files.isdir]).name};
dirnames = dirnames(~ismember(dirnames,{'.', '..'} ));

%First we will count the number of .mat files to create the final matrix
tot = 0;
for i = 1:size(dirnames,2)
    %First count the number of *mat files
    thisdir = char(dirnames(:,i));
    subdirinfo = dir(fullfile(thisdir, '*.mat'));
    tot        = tot + size(subdirinfo, 1);
end

%Next we will open all files and add the vertices to the final matrix
landmark_matrix = zeros(tot, 21480);
landmark_ids    = cell(tot,1);
count = 1;
for i = 1:size(dirnames,2)
    thisdir    = char(dirnames(:,i));
    subdirinfo = dir(fullfile(thisdir, '*.mat'));
    for r = 1:size(subdirinfo,1)
        load(fullfile(subdirinfo(1).folder, subdirinfo(r).name))
        vert = reshape(in.Vertices', [1 21480]);
        landmark_ids(count,1)     = cellstr(erase(subdirinfo(r).name, ".mat"));
        landmark_matrix(count, :) = vert;
        count = count + 1;
        strcat(int2str(count), '.... ', subdirinfo(r).name)
    end
end

landmark_facets = in.Faces;

%Next we will export the landmark matrix
cd(folder.save);
save('landmark_facets', 'landmark_facets')
save('landmark_matrix', 'landmark_matrix')
save('landmark_ids', 'landmark_ids')