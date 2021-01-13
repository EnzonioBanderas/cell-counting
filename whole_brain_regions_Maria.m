
%run SHARPtrack analyze_ROIs first for the variable roi_table
%find groups of strings for sub areas of the brain

% file_name = 'amp';
% [TEXT, everything] = xlsread(file_name);
% roi_table = everything(:,1);
% roi_table = cell2table(roi_table); 
% clear everything
% clear TEXT
processed_folder = uigetdir('', 'Select processed folder');
load(fullfile(processed_folder, 'roi_table_all'))
load(fullfile(processed_folder, 'pixel_table_all'))
load(fullfile(processed_folder, 'pername_table'))

roi_table = roi_table_all ; %change this if necessary 
c = roi_table(:,1); % creates a variable with a list of every cell in a brain region
[d, id] = findgroups(c); % d allocates a group number to every variable in c, ID is the unique name of every group in alphabetical order
d2 = sort(d);
[group_number,ia,ic] = unique(d2);
d_counts = accumarray(ic,1);
sort_d_counts = [group_number, d_counts];
sort_d_counts = sortrows(sort_d_counts,2);
list= string(id.name); % transorms X to string variable
table_brain_areas = table(d_counts, group_number, list);
sorted_table_brain_areas = sortrows(table_brain_areas);
list2 = sorted_table_brain_areas.list;

%%
Y = [sorted_table_brain_areas.d_counts, sorted_table_brain_areas.list]; % concatenates the two variables and shows how many clicked points there were in every brain area

X = categorical(sorted_table_brain_areas.list); % X = categorical(id.name) creates a categorical array from the array id.name. The categories of x are the sorted unique values from Id.name.
X = reordercats(X,sorted_table_brain_areas.list); % reorders the categories based on the name

%% Plotting horizontal bar graph with brain regions 
%plot all the areas
figure
barh(X,sorted_table_brain_areas.d_counts,'FaceColor',[0.82 0.33 0.33],'EdgeColor', 'none', 'BarWidth', .6);
yticklabels(list2)
ylabel('Brain areas')
set(gca,'FontSize', 10, 'Fontname', 'Calibri')
xlabel('C-fos positive neurons')
ylabel('Brain areas')
set(gca,'FontSize', 10, 'Fontname', 'Calibri')
xlabel('C-fos positive neurons')

ax = gca;
co = ax.XAxis.Color;
ax.XAxis.Color = 'black';
ax = gca;
co = ax.YAxis.Color;
ax.YAxis.Color = 'black';

set(gcf,'Color','w')
set(gca,'Color','w')

%title('C-fos positive neurons in the Midbrain, F3_JI')

%legend('Total number of c-fos positive neurons = 14611')
%legend('Location','southoutside','FontSize',10,'TextColor','black'); 
%legend boxoff


%% plot the areas with more neurons (selection)
% index = 98:108; % change here depending on how long the list is, how many areas
nTop = 10;

X_max = X(end-nTop:end,:); 
sorted_table_max = sorted_table_brain_areas(end-nTop:end,:);
list3 = list2(end-nTop:end,:); 
X_max = cellstr(X_max);
X_max = cell2table(X_max);
X_max = table2array(X_max);
X_max = categorical (X_max);

figure
X_max_reordered = reordercats(X_max, cellstr(X_max));
barh(X_max_reordered, sorted_table_max.d_counts,'FaceColor',[0.82 0.33 0.33],'EdgeColor', 'none');
% yticklabels(list3)
ylabel('Brain areas')
xlabel('Neurons')
set(gca,'TickDir','out');
box('off');

%% Select and Plot other regions 

%find groups of strings for bigger areas of the brain
name ={'Substantia nigra'};
k = strfind(Y(:,2),name);
emptytozero = cellfun('isempty',k); 
k(emptytozero) = {0} ;
k = cell2mat(k);

% amount of roi's in bigger brain area
for i=1:length(k)
    if k(i)>=1
       k(i)=1;
    end
end
Y2 = char(Y(:,1));
Y2 = str2num(Y2);
k2=k.'*Y2;
Total=[name,k2];

%%
%find groups of strings for bigger areas of the brain
name ={'Periaqueductal gray'};
k = strfind(Y(:,2),name);
emptytozero = cellfun('isempty',k); 
k(emptytozero) = {0} ;
k = cell2mat(k);

% amount of roi's in bigger brain area
for i=1:length(k)
    if k(i)>=1
       k(i)=1;
    end
end

k2=k.'*Y2;
Total=[Total; name,k2];

%% Plotting selection of areas

Total2 = cell2mat(Total(:,2));
figure
barh(Total2,'FaceColor',[0.82 0.33 0.33],'EdgeColor', 'none', 'BarWidth', .6)
yticklabels(Total(:,1))
set(gca,'FontSize', 12, 'Fontname', 'Calibri')
ylabel('Brain areas')
xlabel('Neurons')
set(gca,'TickDir','out');
box('off');


%Change xlim if necessary!!!!
% xlim([0 400]);


