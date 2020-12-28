% ------------------------------------------------------------------------
%        Get ROI reference-space locations and region annotations
% ------------------------------------------------------------------------


%% SET FILE LOCATIONS OF TRANSFORM AND ROIS

% IMP!! 
% 1- Change the slide number! For proper saving later on.
% 2- Change destination excel file and range.
% 3- Change location of transform and roi array


slicename = 'amp_s23'; 

% file location of transform and transformed image
transform_location = 'D:\Fabian_stainings\20201126_Pax5 staining with signal amplification\pax5_IF_amp\test1\transformations\downsampled_pax5_IF_amp_s23c1-3_transform_data';
transformed_slice_location = 'D:\Fabian_stainings\20201126_Pax5 staining with signal amplification\pax5_IF_amp\test1\transformations\downsampled_pax5_IF_amp_s23c1-3_transformed.tif';

% Using a set of x and y coordinates from the ImageJ function Analyze
% Particles to generate an ROI image
roi_array = zeros(800,1140,'uint8');
roi_array_values = csvread('D:\Fabian_stainings\20201126_Pax5 staining with signal amplification\pax5_IF_amp\test1\transformations\downsampled_pax5_IF_amp_s23c1-3_transformed.csv',1,1);
x = roi_array_values(:, 2); % IMPORTANT: the x and y coordinates have to be changed bc the images are 1140x800 and the reference atlas is 800x1140.
y = roi_array_values(:, 1); %is the ix and y are not changed, everything is flipped around. 

for i = 1:length(roi_array_values)-1
  roi_array(x(i),y(i)) = 200; 
end
    
% directory of reference atlas files
annotation_volume_location = 'D:\cfos\SharpTrack\annotation_volume_10um_by_index.npy';
structure_tree_location = 'D:\cfos\SharpTrack\structure_tree_safe_2017.csv';

%% LOAD THE DATA

% load the transformed slice image
transformed_slice_image = imread(transformed_slice_location);

% load the transform from the transform file
transform_data = load(transform_location);
transform_data = transform_data.save_transform;

% get the actual transformation from slice to atlas
slice_to_atlas_transform = transform_data.transform;

% get the position within the atlas data of the transformed slice
slice_num = transform_data.allen_location{1};
slice_angle = transform_data.allen_location{2};

% load the rois rois = imread(roi_location);

% if the rois come from a transformed roi image of non-contiguous roi
% pixels (e.g. an ROI pixel for each neuron), then run this line to ensure
% a one-to-one mapping between ROIs in the original and transformed images:
%rois = uint8(imregionalmax(rois));
roi_array = uint8(imregionalmax(roi_array));

% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end


%% GET REFERENCE-SPACE LOCATIONS AND REGION ANNOTATIONS FOR EACH ROI

% I will do this for every *pixel* in the roi image with nonzero value, but
% this code can be modified, e.g. to do it by clusters of pixels


% show the transformed ROI, together with the transformed image
%figure; imshow(imfuse(rois, transformed_slice_image));

figure; 
imshow(imfuse(roi_array, transformed_slice_image));
title('transformed slice image, fused with ROIs')

figure;
imshow(imfuse(roi_array, roi_array))

% make sure the rois are in a properly size image
%assert(size(rois,1)==800&size(rois,2)==1140&size(rois,3)==1,'roi image is
%not the right size');
assert(size(roi_array,1)==800&size(roi_array,2)==1140&size(roi_array,3)==1,'roi image is not the right size');

% initialize array of locations (AP, DV, ML relative to bregma) in
% reference space and the correponding region annotations
%roi_location = zeros(sum(rois(:)>0),3); roi_annotation =
%cell(sum(rois(:)>0),3);
roi_location = zeros(sum(roi_array(:)>0),3);
roi_annotation = cell(sum(roi_array(:)>0),3);

% get location and annotation for every roi pixel
[pixels_row, pixels_column] = find(roi_array>0);

% generate other necessary values
bregma = allenCCFbregma(); % bregma position in reference data space
atlas_resolution = 0.010; % mm
offset_map = get_offset_map(slice_angle);

% loop through every pixel to get ROI locations and region annotations
for pixel = 1:length(pixels_row)
    
    % get the offset from the AP value at the centre of the slice, due to
    % off-from-coronal angling
    offset = offset_map(pixels_row(pixel),pixels_column(pixel));
    
    % use this and the slice number to get the AP, DV, and ML coordinates
    ap = -(slice_num-bregma(1)+offset)*atlas_resolution;
    dv = (pixels_row(pixel)-bregma(2))*atlas_resolution;
    ml = (pixels_column(pixel)-bregma(3))*atlas_resolution;

    roi_location(pixel,:) = [ap dv ml];
    
    % finally, find the annotation, name, and acronym of the current ROI
    % pixel
    ann = av(slice_num+offset,pixels_row(pixel),pixels_column(pixel));
    name = st.safe_name{ann};
    acr = st.acronym{ann};
    
    roi_annotation{pixel,1} = ann;
    roi_annotation{pixel,2} = name;
    roi_annotation{pixel,3} = acr;

end
 
 roi_table = table(roi_annotation(:,2),roi_annotation(:,3), ...
                    roi_location(:,1),roi_location(:,2),roi_location(:,3), roi_annotation(:,1), ...
     'VariableNames', {'name', 'acronym', 'AP_location', 'DV_location', 'ML_location', 'avIndex'});
 
%% Save roi_table as a .mat file and write roi_table information in an excel file. 

save(slicename, 'roi_table'); 
%writetable(roi_table,'amp.xlsx');
%writetable(roi_table, fullfile(pwd, 'amp.xlsx'),'Range','A2374');



