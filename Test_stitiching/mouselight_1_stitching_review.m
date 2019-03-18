%% Path setup 
raw_data_grid = load('mouselight_1_raw_data_info.mat');
pipelineoutputfolder = '/nrs/mouselight/pipeline_output/2018-08-15_pipeline_test';
line_fix_folder = fullfile(pipelineoutputfolder, 'stage_1_line_fix_output');
descriptor_folder = fullfile(pipelineoutputfolder, 'stage_2_descriptor_output');
point_match_folder = fullfile(pipelineoutputfolder, 'stage_3_point_match_output');
%% Debug stitching
test_grid_xyz = [232 28 1540];
test_grid_idx = find(all(raw_data_grid.stage_grid_xyz == test_grid_xyz, 2));
test_tile_fp = raw_data_grid.image_filepath{test_grid_idx};
% test_image = deployedtiffread(test_tile_fp);
% implay(test_image); 
%%
% Load descriptor and matched descriptor
clearvars test_desc match_x match_y match_z
dir_part = strsplit(test_tile_fp, '/');
relative_path = dir_part(end-3:end-1);
relative_path = fullfile(relative_path{:});
descriptor_fp = dir(fullfile(descriptor_folder, relative_path, '*-desc.0.mat'));
descriptor_fp = fullfile(descriptor_fp.folder, descriptor_fp.name);
figure('Visible', 'on')
clf

if ~isempty(descriptor_fp) && isfile(descriptor_fp)
    test_desc = load(descriptor_fp);
end
point_match_x_fp = dir(fullfile(point_match_folder, relative_path, '*X.mat'));
if ~isempty(point_match_x_fp) 
    point_match_x_fp = fullfile(point_match_x_fp.folder, point_match_x_fp.name);
    if isfile(point_match_x_fp)
        match_x = load(point_match_x_fp);
    else
        disp('X-matching file does not exist');
    end
    if ~isempty(match_x.tmp_xy_match_result.paireddescriptor.X)
        subplot(2,2,2);
        [sub_1, sub_2] = deal([]);
        sub_1 = match_x.tmp_xy_match_result.paireddescriptor.X;
        sub_2 = match_x.tmp_xy_match_result.paireddescriptor.Y;
        scatter3(sub_1(:,2), sub_1(:,1), sub_1(:,3));
        hold on
        scatter3(sub_2(:,2), sub_2(:,1), sub_2(:,3));
        title('Match in X direction');
        daspect([1,1,1]);
    else
        disp('X-matching is empty');
    end
else
    disp('X_matching file does not exist');
end

point_match_y_fp = dir(fullfile(point_match_folder, relative_path, '*Y.mat'));
if ~isempty(point_match_y_fp) 
    point_match_y_fp = fullfile(point_match_y_fp.folder, point_match_y_fp.name);
    if isfile(point_match_y_fp)
        match_y = load(point_match_y_fp);
    else
        disp('Y_matching file does not exist');
    end
    if ~isempty(match_y.tmp_xy_match_result.paireddescriptor.X)
        subplot(2,2,3);
        [sub_1, sub_2] = deal([]);
        sub_1 = match_y.tmp_xy_match_result.paireddescriptor.X;
        sub_2 = match_y.tmp_xy_match_result.paireddescriptor.Y;
        scatter3(sub_1(:,2), sub_1(:,1), sub_1(:,3));
        hold on
        scatter3(sub_2(:,2), sub_2(:,1), sub_2(:,3));
        title('Match in Y direction');
        daspect([1,1,1]);
    else
        disp('Y-matching is empty');
    end
else
    disp('Y_matching file does not exist');
end
point_match_z_fp = dir(fullfile(point_match_folder, relative_path, '*Z.mat'));
point_match_z_1_fp = dir(fullfile(point_match_folder, relative_path, '*Z-1.mat'));
if ~isempty(point_match_z_fp) 
    point_match_z_fp = fullfile(point_match_z_fp.folder, point_match_z_fp.name);
    point_match_z_1_fp = fullfile(point_match_z_1_fp.folder, point_match_z_1_fp.name);
    if isfile(point_match_z_fp)
        match_z = load(point_match_z_fp);
        match_z_1 = load(point_match_z_1_fp);
        disp(match_z.paireddescriptor);
        disp(match_z_1.paireddescriptor);
    else
        disp('Z-matching file does not exist');
    end
    if ~isempty(match_z.paireddescriptor.X)
        subplot(2,2,1)
        [sub_1, sub_2] = deal([]);
        sub_1 = match_z.paireddescriptor.X;
        sub_2 = match_z.paireddescriptor.Y;
        scatter3(sub_1(:,2), sub_1(:,1), sub_1(:,3));
        hold on
        scatter3(sub_2(:,2), sub_2(:,1), sub_2(:,3));
        title('Match in Z direction');
        daspect([1,1,1]);
    else
        disp('Z-matching is empty');
    end
else
    disp('Z_matching file does not exist');
end
% Visualization 
subplot(2,2,4);
[sub_1, sub_2] = deal([]);
sub_1 = test_desc.skl_sub;
sub_2 = test_desc.edge_sub;
if ~isempty(sub_1)
    scatter3(sub_1(:,2), sub_1(:,1), sub_1(:,3));
    hold on 
else
    disp('Skeleton is empty');
end
if ~isempty(sub_2)
    if size(sub_2, 1) > 100000
        scatter3(sub_2(1:10:end,2), sub_2(1:10:end,1), sub_2(1:10:end,3));
    else
        scatter3(sub_2(:,2), sub_2(:,1), sub_2(:,3));
    end
else
    disp('Edge is empty');
end
daspect([1,1,1]);
% New z-matching
figure;
scatter3(match_z_1.paireddescriptor.X(:,2), match_z_1.paireddescriptor.X(:,1), ...
    match_z_1.paireddescriptor.X(:,3));
hold on 
scatter3(match_z_1.paireddescriptor.Y(:,2), match_z_1.paireddescriptor.Y(:,1), ...
    match_z_1.paireddescriptor.Y(:,3));
%% Update feautres with lower global threshold 

%% Check if the skeleton is from the image
% Load image
test_image = deployedtiffread(test_tile_fp);
test_image_proj = max(test_image, [], 3);
% Load new descriptor
skl_sub_2d = test_desc.skl_sub(:, 1:2);
figure;
scatter(skl_sub_2d(:,2), skl_sub_2d(:,1));