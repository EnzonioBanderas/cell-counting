% ------------------------------------------------------------------------
%        Get ROI reference-space locations and region annotations from cell
%        locations in csv file (Image J)
% ------------------------------------------------------------------------

clearvars; 
% Specify the folder where the files live 
% Folder = 'D:\Fabian_stainings\20201126_Pax5 staining with signal amplification\pax5_IF_amp\test2\transformations'; 
transformations_folder = uigetdir('', 'Select transformations folder');
processed_folder = fullfile(transformations_folder, '..');
roifused_folder = fullfile(processed_folder, 'roifused_folder');
mkdir(roifused_folder)

% Get a list of all files in the folder with the desired file name pattern and open them. 
filePattern_transfrom = fullfile(transformations_folder, '*.mat'); 
filePattern_slice = fullfile(transformations_folder, '*.tif'); 
filePattern_roi = fullfile(transformations_folder, '*.csv');

Files_transfrom = dir(filePattern_transfrom);
Files_slice = dir(filePattern_slice);
Files_roi = dir(filePattern_roi); 

roi_table = {}; %preallocation of the array that will be filled with every loop iteration (info from every slice)

pixel_table = cell(length(Files_transfrom), 1);
roi_table = cell(length(Files_transfrom), 1);
for k = 1 : length(Files_transfrom)
    % file location of transform, transformed image and csv file fro ROI creation
    baseFileName_transform = Files_transfrom(k).name;
    transform_location = fullfile(Files_transfrom(k).folder, baseFileName_transform);
    fprintf(1, 'Now reading %s\n', baseFileName_transform);

    baseFileName_slice = Files_slice(k).name; 
    transformed_slice_location = fullfile(Files_slice(k).folder, baseFileName_slice);

    baseFileName_slice = Files_roi(k).name; 
    roi_array_values = csvread((fullfile(Files_slice(k).folder, baseFileName_slice)),1,1);
    % Using a set of x and y coordinates from the ImageJ function Analyze Particles to generate an ROI image
    roi_array = zeros(800,1140,'uint8');
    x = roi_array_values(:, 2); % IMPORTANT: the x and y coordinates have to be changed bc the images are 1140x800 and the reference atlas is 800x1140.
    y = roi_array_values(:, 1); %is the ix and y are not changed, everything is flipped around. 

    for i = 1:length(roi_array_values)-1
        roi_array(x(i),y(i)) = 200; 
    end

    % directory of reference atlas files
    % annotation_volume_location = 'D:\cfos\SharpTrack\annotation_volume_10um_by_index.npy';
    % structure_tree_location = 'D:\cfos\SharpTrack\structure_tree_safe_2017.csv';
    annotation_volume_location = fullfile('SharpTrack', 'annotation_volume_10um_by_index.npy');
    structure_tree_location = fullfile('SharpTrack', 'structure_tree_safe_2017.csv');

    % load the transformed slice image
    transformed_slice_image = imread(transformed_slice_location);
    transformed_slice_image = max(transformed_slice_image, [], 3); % take the maximum over RGB dim, assume only one channel contains all info

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


    % GET REFERENCE-SPACE LOCATIONS AND REGION ANNOTATIONS FOR EACH ROI

    % Do this for every *pixel* in the roi image with nonzero value, but
    % this code can be modified, e.g. to do it by clusters of pixels

    FIG = figure; 
    imshow(imfuse(roi_array, transformed_slice_image));
    title('transformed slice image, fused with ROIs')
    saveas(FIG, fullfile(roifused_folder, ['roifused_', 's', sprintf('%02d', k)]), 'tif')

%     figure;
%     imshow(imfuse(roi_array, roi_array))

    % make sure the rois are in a properly size image
    %assert(size(rois,1)==800&size(rois,2)==1140&size(rois,3)==1,'roi image is
    %not the right size');
    assert(size(roi_array,1)==800&size(roi_array,2)==1140&size(roi_array,3)==1,'roi image is not the right size');

    % initialize array of locations (AP, DV, ML relative to bregma) in
    % reference space and the correponding region annotations
    roi_location = zeros(sum(roi_array(:)>0),3);
    roi_annotation = cell(sum(roi_array(:)>0),3);

    % get location and annotation for every roi pixel
    [pixels_row, pixels_column] = find(roi_array>0);

    % generate other necessary values
    bregma = allenCCFbregma(); % bregma position in reference data space
    atlas_resolution = 0.010; % mm
    offset_map = get_offset_map(slice_angle);

    % loop through every pixel in all images and count number of pixels and summed
    % signal intensity per area
    % NEED: signal arrays for each image
    pixel_location = zeros(numel(transformed_slice_image), 3);
    signal_vec = zeros(numel(transformed_slice_image), 1);
    pixel_annotation = cell(numel(transformed_slice_image), 3);
    iPixel = 0;
    for iX = 1:size(transformed_slice_image, 1)
        for iY = 1:size(transformed_slice_image, 2)
            iPixel = iPixel + 1;
            % Create a table of names (and acronyms and annotations) and
            % signal values. Calculate mean signal intensity per ROI in
            % addition to number of pixels. Use number of pixels to
            % normalize results.
            signal_vec(iPixel) = transformed_slice_image(iX, iY);

            % get the offset from the AP value at the centre of the slice, due to
            % off-from-coronal angling
            offset = offset_map(iX, iY);

            % use this and the slice number to get the AP, DV, and ML coordinates
            ap = -(slice_num-bregma(1)+offset)*atlas_resolution;
            dv = (iX-bregma(2))*atlas_resolution;
            ml = (iY-bregma(3))*atlas_resolution;

            pixel_location(iPixel,:) = [ap dv ml];

            % finally, find the annotation, name, and acronym of the current ROI
            % pixel
            ann = av(slice_num+offset, iX, iY);
            name = st.safe_name{ann};
            acr = st.acronym{ann};

            pixel_annotation{iPixel,1} = ann;
            pixel_annotation{iPixel,2} = name;
            pixel_annotation{iPixel,3} = acr;
        end
    end
    all_pixel_table = table(pixel_annotation(:,2), pixel_annotation(:,3), ...
                pixel_location(:,1), pixel_location(:,2), pixel_location(:,3), pixel_annotation(:,1), ...
                signal_vec, ...
    'VariableNames', {'name', 'acronym', 'AP_location', 'DV_location', 'ML_location', 'avIndex', 'signal'});
    [uniq, uniq_idx, or_idx] = unique(all_pixel_table.name);
    pixel_count = accumarray(or_idx(:), 1, [], @sum);
    signal_sum = accumarray(or_idx(:), all_pixel_table.signal, [], @sum);
    
    pixel_table{k} = table(uniq, all_pixel_table.acronym(uniq_idx), ...
                all_pixel_table.avIndex(uniq_idx), ...
                pixel_count, signal_sum, ...
    'VariableNames', {'name', 'acronym', 'avIndex', 'pixel_count', 'signal_sum'});

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
        ann = av(slice_num+offset,pixels_row(pixel), pixels_column(pixel));
        name = st.safe_name{ann};
        acr = st.acronym{ann};

        roi_annotation{pixel,1} = ann;
        roi_annotation{pixel,2} = name;
        roi_annotation{pixel,3} = acr;

    end

    roi_table{k} = table(roi_annotation(:,2),roi_annotation(:,3), ...
                roi_location(:,1),roi_location(:,2),roi_location(:,3), roi_annotation(:,1), ...
    'VariableNames', {'name', 'acronym', 'AP_location', 'DV_location', 'ML_location', 'avIndex'});

    close(FIG);

end

roi_table_all = vertcat(roi_table{:});
save(fullfile(processed_folder, 'roi_table_all'), 'roi_table_all'); %change first input to whatever name you want 
pixel_table_all = vertcat(pixel_table{:});
save(fullfile(processed_folder, 'pixel_table_all'), 'pixel_table_all'); %change first input to whatever name you want 

% count number of ROIs per name and number of total pixels and total
% intensity (probably per mouse normalization necessary)
[uniq, uniq_idx, or_idx] = unique(roi_table_all.name);
roi_count = accumarray(or_idx(:), 1, [], @sum);
pername_roi_table = table(uniq, roi_table_all.acronym(uniq_idx), ...
            roi_table_all.avIndex(uniq_idx), ... 
            roi_count, ...
            'VariableNames', {'name', 'acronym', 'avIndex', ...
            'roi_count'});
[uniq, uniq_idx, or_idx] = unique(pixel_table_all.name);
signal_sum = accumarray(or_idx(:), pixel_table_all.signal_sum, [], @sum);
pixel_count = accumarray(or_idx(:), pixel_table_all.pixel_count, [], @sum);
signal_mean = signal_sum ./ pixel_count;
pername_pixel_table = table(uniq, pixel_table_all.acronym(uniq_idx), ...
            pixel_table_all.avIndex(uniq_idx), ... 
            signal_sum, pixel_count, signal_mean, ...
            'VariableNames', {'name', 'acronym', 'avIndex', ...
            'signal_sum', 'pixel_count', 'signal_mean'});

% merge collapsed ROI table with pixel_table to create pername table
pername_table = outerjoin(pername_pixel_table, pername_roi_table, 'Keys', 'name', 'RightVariables', 'roi_count');
pername_table.roi_count(isnan(pername_table.roi_count)) = 0;
pername_table.roi_count_perpixel = pername_table.roi_count ./ pername_table.pixel_count;
[~, sort_idx] = sort(pername_table.roi_count_perpixel, 'descend');
pername_table = pername_table(sort_idx, :);
save(fullfile(processed_folder, 'pername_table'), 'pername_table');

% function location2meta(xPixel, yPixel, offset_map, slice_num)
%     % get the offset from the AP value at the centre of the slice, due to
%     % off-from-coronal angling
%     offset = offset_map(xPixel, yPixel);
% 
%     % use this and the slice number to get the AP, DV, and ML coordinates
%     ap = -(slice_num-bregma(1)+offset)*atlas_resolution;
%     dv = (pixels_row(pixel)-bregma(2))*atlas_resolution;
%     ml = (pixels_column(pixel)-bregma(3))*atlas_resolution;
% 
%     roi_location(pixel,:) = [ap dv ml];
% 
%     % finally, find the annotation, name, and acronym of the current ROI
%     % pixel
%     ann = av(slice_num+offset,pixels_row(pixel), pixels_column(pixel));
%     name = st.safe_name{ann};
%     acr = st.acronym{ann};
% end



