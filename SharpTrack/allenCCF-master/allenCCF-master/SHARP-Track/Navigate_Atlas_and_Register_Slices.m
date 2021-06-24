% ------------------------------------------------------------------------
%          Run Allen Atlas Browser
% ------------------------------------------------------------------------


%% ENTER FILE LOCATION AND PROBE-SAVE-NAME


% directory of histology

% processed_images_folder = 'D:\Fabian_stainings\Exported pictures of WT brain for Maria\2'; 
processed_images_folder = uigetdir('' , 'select processed images folder');

% name the saved probe points, to avoid overwriting another set of probes going in the same folder
probe_save_name_suffix = ''; 

% directory of reference atlas filesa
annotation_volume_location = fullfile('SharpTrack', 'annotation_volume_10um_by_index.npy');
structure_tree_location = fullfile('SharpTrack', 'structure_tree_safe_2017.csv');
template_volume_location = fullfile('SharpTrack', 'template_volume_10um.npy');

% plane to view ('coronal', 'sagittal', 'transverse')
plane = 'coronal';

%% GET PROBE TRAJECTORY POINTS

% load the reference brain and region annotations
if ~exist('av','var') || ~exist('st','var') || ~exist('tv','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
    tv = readNPY(template_volume_location);
end

% select the plane for the viewer
if strcmp(plane,'coronal')
    av_plot = av;
    tv_plot = tv;
elseif strcmp(plane,'sagittal')
    av_plot = permute(av,[3 2 1]);
    tv_plot = permute(tv,[3 2 1]);
elseif strcmp(plane,'transverse')
    av_plot = permute(av,[2 3 1]);
    tv_plot = permute(tv,[2 3 1]);
end

% create Atlas viewer figure
f = figure('Name','Atlas Viewer'); 

% show histology in Slice Viewer
try; figure(slice_figure_browser); title('');
catch; slice_figure_browser = figure('Name','Slice Viewer'); end
reference_size = size(tv_plot);
sliceBrowser(slice_figure_browser, processed_images_folder, f, reference_size);


% use application in Atlas Transform Viewer
% use this function if you have a processed_images_folder with appropriately processed .tif histology images
[f, ud] = AtlasTransformBrowser(f, tv_plot, av_plot, st, slice_figure_browser, processed_images_folder, probe_save_name_suffix); 


% use the simpler version, which does not interface with processed slice images
% just run these two lines instead of the previous 5 lines of code
% 
% save_location = processed_images_folder;
% f = allenAtlasBrowser(tv_plot, av_plot, st, save_location, probe_save_name_suffix);
% 

figure
imshow(curr_slice_trans);
hold on

% OverlayImage=[];
% F = scatteredInterpolant(Y, X, curr_slice_trans,'linear');
% for i = 1:height-1
%    for j = 1:width-1
%           OverlayImage(i,j) = F(i,j);
%    end
% end
% alpha = (~isnan(OverlayImage))*0.6;

OverlayImage = imshow( ref );
caxis auto
colormap( OverlayImage.Parent, gray );
colorbar( OverlayImage.Parent );
curr_size = size(ref);
set( OverlayImage, 'AlphaData', ones(curr_size(1:2))*0.6 );




% Enter data by hand if data from a ThingSpeak channel is not available.
strength = [-90 -90 -90 -90 -40 -20 -22.4 -45 -35 -41 -44 -55 -40 -75 -26]';
% Read data from a ThingSpeak channel.
% Uncomment the next line to read from ThingSpeak.
% strength = thingSpeakRead(CHANNEL_ID, 'ReadKey',READ_API_KEY,'numPoints',15,'fields',FIELD_NUMBER');
X = [10 550 550 10 50 234 393 129 237 328 448 225 344 457 477]';
Y = [10 10 410 410 293 210 202 132 130 142 141 272 268 274 200]';

strengthPercent = 2*(strength+100)/100;

picture = imread('https://www.mathworks.com/help/examples/thingspeak/win64/CreateHeatmapOverlayImageTSExample_02.png'); 
[height,width,depth] = size(picture);

OverlayImage=[];
F = scatteredInterpolant(Y, X, strengthPercent,'linear');
for i = 1:height-1
   for j = 1:width-1
          OverlayImage(i,j) = F(i,j);
   end
end
alpha = (~isnan(OverlayImage))*0.6;

imshow(picture);
hold on
OverlayImage = imshow( OverlayImage );
caxis auto  
colormap( OverlayImage.Parent, jet );
colorbar( OverlayImage.Parent );
set( OverlayImage, 'AlphaData', alpha );
