% ------------------------------------------------------------------------
%        Get ROI reference-space locations and region annotations from cell
%        locations in csv file (Image J)
% ------------------------------------------------------------------------

clearvars; 
% Specify the folder where the files live 
% Folder = 'D:\Fabian_stainings\20201126_Pax5 staining with signal amplification\pax5_IF_amp\test2\transformations'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transformations_folder = uigetdir('', 'Select transformations folder');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_path = '/mnt/tosh/Projects/MEP/cell-counting/Immunostainings_Pax5_mergedSelectionMariaSS';
tfa_struct = dir(data_path);
tfa_struct = tfa_struct([tfa_struct.isdir]);
nTfa = length(tfa_struct);
tfa = cell(1, nTfa);
for iTfa = 3:nTfa
    tfa{iTfa} = fullfile(data_path, tfa_struct(iTfa).name, '/processed/transformations');
end
tfa= tfa(3:end);
% tfa = {'/mnt/tosh/Projects/MEP/cell-counting/Immunostainings_Pax5_extra/20210320_Pax5++_2A#07_M3_DAPI,TH-AF488,Pax5-Ampl-AF594-01/processed/transformations', ...
%     '/mnt/tosh/Projects/MEP/cell-counting/Immunostainings_Pax5_extra/20210320_Pax5+-_2A#02_M1_DAPI,TH-AF488,Pax5-Ampl-AF594-01/processed/transformations', ...
%     '/mnt/tosh/Projects/MEP/cell-counting/Immunostainings_Pax5_extra/20210320_Pax5+-_2A#05_M2_DAPI,TH-AF488,Pax5-Ampl-AF594-01/processed/transformations', ...
%     '/mnt/tosh/Projects/MEP/cell-counting/Immunostainings_Pax5_extra/20210322_Pax5-pR31Q-_2A#09_M4_DAPI,TH-AF488,Pax5-Ampl-AF594-01/processed/transformations', ...
%     '/mnt/tosh/Projects/MEP/cell-counting/Immunostainings_Pax5_extra/20210322_Pax5-pR31Q-_2B#01_M5_DAPI,TH-AF488,Pax5-Ampl-AF594-01/processed/transformations', ...
%     '/mnt/tosh/Projects/MEP/cell-counting/Immunostainings_Pax5_extra/20210322_Pax5-pR31Q-_2B#03_M6_DAPI,TH-AF488,Pax5-Ampl-AF594-01/processed/transformations', ...
%     '/mnt/tosh/Projects/MEP/cell-counting/Immunostainings_Pax5_extra/20210323_Pax5+-_2B#06_M7_DAPI,TH-AF488,Pax5-Ampl-AF594-01/processed/transformations', ...
%     '/mnt/tosh/Projects/MEP/cell-counting/Immunostainings_Pax5_extra/20210325_Pax5++_2B#08_M8_DAPI,TH-AF488,Pax5-Ampl-AF594-01/processed/transformations', ...
%     '/mnt/tosh/Projects/MEP/cell-counting/Immunostainings_Pax5_extra/20210325_Pax5++_2B#10_M9_DAPI,TH-AF488,Pax5-Ampl-AF594-01/processed/transformations'};
for iTfa = 1:length(tfa)

transformations_folder = tfa{iTfa};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

processed_folder = fullfile(transformations_folder, '..');

% Get a list of all files in the folder with the desired file name pattern and open them. 
filePattern_slice = fullfile(transformations_folder, '*.tif'); 
filePattern_transfrom = fullfile(transformations_folder, '*.mat'); 
Files_slice = dir(filePattern_slice);
Files_transfrom = dir(filePattern_transfrom);

% Load reference
annotation_volume_location = fullfile('SharpTrack', 'annotation_volume_10um_by_index.npy');
structure_tree_location = fullfile('SharpTrack', 'structure_tree_safe_2017.csv');
disp('loading reference atlas...')
av = readNPY(annotation_volume_location);
st = loadStructureTree(structure_tree_location);

for k = 1 : length(Files_slice)
    % file location of transform, transformed image and csv file fro ROI creation
    baseFileName_transform = Files_transfrom(k).name;
    transform_location = fullfile(Files_transfrom(k).folder, baseFileName_transform);
    fprintf(1, 'Now reading %s\n', baseFileName_transform);
    baseFileName_slice = Files_slice(k).name; 
    transformed_slice_location = fullfile(Files_slice(k).folder, baseFileName_slice);
    
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
%     transformed_slice_annotation = transformed_slice_image;
%     transformed_slice_annotation = zeros(size(transformed_slice_image));
    transformed_slice_annotation = zeros(size(transformed_slice_image), 'int16');
    for iX = 1:size(transformed_slice_image, 1)
%         fprintf([num2str(iX), ' \n'])
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
%             name = st.safe_name{ann};
%             acr = st.acronym{ann};

            pixel_annotation{iPixel,1} = ann;
            
%             pixel_annotation{iPixel,2} = name;
%             pixel_annotation{iPixel,3} = acr;
%             
%             pixel_annotation{iPixel, 2} = iX;
%             pixel_annotation{iPixel, 3} = iY;
            
            % assign annotation to slice
            transformed_slice_annotation(iX, iY) = ann;
            
        end
    end
    
    % isolate in matlab as mask
    % loop and load in imagej and multiply and threshold
    % save results and combine post loop
    
%     figure;
%     imshow(transformed_slice_image)    
%     figure;
%     imshow(transformed_slice_annotation)
%     caxis([min([pixel_annotation{:,1}]), max([pixel_annotation{:,1}])])
    
    % Loop through annotation
    uniqAnn = unique(transformed_slice_annotation(:));
    uniqAnn_length = length(uniqAnn);
%     transformed_slice_annotation_stack = zeros([size(transformed_slice_annotation), uniqAnn_length]);
    transformed_slice_stack_folder = fullfile(Files_slice(k).folder, 'structure_stacks');
    mkdir(transformed_slice_stack_folder)
    baseFileName_slice_split = strsplit(baseFileName_slice, '.');
    transformed_slice_stack_location  = fullfile(transformed_slice_stack_folder, [baseFileName_slice_split{1}, '_structureStack.tif']);
    for iAnn = 1:uniqAnn_length
        transformed_slice_annotation_stackSlice = zeros(size(transformed_slice_annotation));
        transformed_slice_annotation_stackSlice(transformed_slice_annotation == uniqAnn(iAnn)) = 1;
%         transformed_slice_annotation_stack(:, :, iAnn) = transformed_slice_annotation_stackSlice;
        
        if iAnn == 1
            imwrite(transformed_slice_annotation_stackSlice, transformed_slice_stack_location)
        else
            imwrite(transformed_slice_annotation_stackSlice, transformed_slice_stack_location, 'WriteMode', 'append')
        end
    end
    
    
    
%     offset = offset_map(iX, iY);
%     % finally, find the annotation, name, and acronym of the current ROI
%     % pixel
%     ann = av(slice_num+offset, iX, iY);
%     name = st.safe_name{ann};
%     acr = st.acronym{ann};
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
