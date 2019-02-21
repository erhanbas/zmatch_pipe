%% To do list
% 1. Determine which region to start with. This region should be: 
%   a. Contains both large surface vessels and capillaries
%   b. Convenient to test and handled locally 
% 2. Run the line shift correction to the selected dataset
% 3. Run skelDescriptor.m in the pipeline-descriptor-master
% 4. Run the featmach
% 5. Run the stitching-master
% 6. Examne the region of bad matching, the data structure of the stitching
% output 
% 7. Improve the stitching
%
%
%
% Read all the raw data position information 
raw_data_root = '/nfs/birdstore-brainbucket2/Vessel/WholeBrain/mouselight_1/Raw_Green';
% Gather raw images file paths
raw_data_info = dir(fullfile(raw_data_root, '**/*.0.tif'));
num_file = numel(raw_data_info);
% Test if the fov_x_size_um and x_size_um are redundant 
tic
parfor tile_idx = 1 : numel(raw_data_info)
    tmp_stage_info = fun_io_load_microscope_info(raw_data_info(tile_idx).folder);
    tmp_field_name  = fieldnames(tmp_stage_info);
    for field_idx = 1 : numel(tmp_field_name)
        raw_data_info(tile_idx).(strip(tmp_field_name{field_idx})) = tmp_stage_info.(tmp_field_name{field_idx});
    end
end
toc
%% Determine the grid position of the valid tile
raw_data_grid = struct;
grid_sub_1 = [raw_data_info.y];
grid_sub_2 = [raw_data_info.x];
grid_sub_3 = [raw_data_info.z];
stage_x_um = [raw_data_info.x_mm];
stage_y_um = [raw_data_info.y_mm];
stage_z_um = [raw_data_info.z_mm];
grid_origin = [min(grid_sub_1), min(grid_sub_2), min(grid_sub_3)];
grid_sub_max = [max(grid_sub_1), max(grid_sub_2), max(grid_sub_3)];
grid_size = grid_sub_max - grid_origin + 1;
grid_sub = bsxfun(@minus, [grid_sub_1;grid_sub_2;grid_sub_3], grid_origin') + 1;

grid_imaged_linear_idx = zeros(grid_size);
raw_data_grid.stage_xyz_um = zeros([3, grid_size]);
raw_data_grid.grid_size = grid_size;
raw_data_grid.num_tile = num_file;
raw_data_grid.image_filepath_array = cell(raw_data_grid.grid_size);

for tile_idx = 1 : num_file
    if ~grid_imaged_linear_idx(grid_sub(1,tile_idx), grid_sub(2, tile_idx), grid_sub(3, tile_idx))
        grid_imaged_linear_idx(grid_sub(1,tile_idx), grid_sub(2, tile_idx), grid_sub(3, tile_idx)) = tile_idx;
        raw_data_grid.image_filepath_array{grid_sub(1,tile_idx), grid_sub(2, tile_idx), grid_sub(3, tile_idx)} = fullfile(raw_data_info(tile_idx).folder, ...
            raw_data_info(tile_idx).name);
        raw_data_grid.stage_xyz_um(:, grid_sub(1,tile_idx), grid_sub(2, tile_idx), grid_sub(3, tile_idx)) = [stage_x_um(tile_idx),...
            stage_y_um(tile_idx), stage_z_um(tile_idx)];
    else
        error('Duplicated tile grid coordinate');
    end
end
raw_data_grid.linear_idx_array = grid_imaged_linear_idx;
raw_data_grid.grid_pos_ind = find(raw_data_grid.linear_idx_array>0);
raw_data_grid.grid_pos_sub = fun_ind2sub(grid_size, raw_data_grid.grid_pos_ind);
raw_data_grid.image_filepath = raw_data_grid.image_filepath_array(raw_data_grid.grid_pos_ind);

raw_data_grid.image_filepath_array(raw_data_grid.grid_pos_ind) = raw_data_grid.image_filepath;
raw_data_grid.FOV_size_um = [raw_data_info(1).fov_y_size_um, raw_data_info(1).fov_x_size_um, 251];
raw_data_grid.block_size = [1536, 1024, 251];
raw_data_grid.voxel_size = raw_data_grid.FOV_size_um ./ raw_data_grid.block_size;
raw_data_grid.overlap_size_um = [raw_data_info(1).fov_y_overlap_um, ...
    raw_data_info(1).fov_x_overlap_um, raw_data_info(1).fov_z_overlap_um];
raw_data_grid.overlap_size = round(raw_data_grid.overlap_size_um ./ raw_data_grid.voxel_size);
save('~/Documents/Github/VesselReconstruction/Metadata/mouselight_1_raw_data_info.mat', '-struct', 'raw_data_grid');
%% Determine region of interest for testing stitching algorithm
% Region with both large vessels and capillaries: 
clc;clear;
raw_data_grid = load('~/Documents/Github/VesselReconstruction/Metadata/mouselight_1_raw_data_info.mat');
DataManager = FileManager;
dataset_name = 'ML_stitching';
for iter_tile = 1 : raw_data_grid.num_tile
    if ~isempty(strfind(raw_data_grid.image_filepath{iter_tile}, '2018-08-23/01/01567'))
        tile_idx = iter_tile;
        test_region_grid_coordinate = raw_data_grid.grid_pos_sub(iter_tile, :);
    end
end

for iter_tile = 1 : prod(raw_data_grid.grid_size)
    if ~isempty(strfind(raw_data_grid.image_filepath_array{iter_tile}, '2018-08-23/01/01567'))
        tile_ind = iter_tile;
%         test_region_grid_coordinate = raw_data_grid.grid_pos_sub(:, iter_tile)';
    end
end

test_region_grid_coordinate = [6, 19, 58];
stack = [num2str(test_region_grid_coordinate, '%02d_'), 'cube'];
num_neighbor = 3;
test_region_file_list = raw_data_grid.image_filepath_array((test_region_grid_coordinate(1)-num_neighbor):(test_region_grid_coordinate(1)+num_neighbor),...
    (test_region_grid_coordinate(2)-num_neighbor):(test_region_grid_coordinate(2)+num_neighbor), ...
    (test_region_grid_coordinate(3)-num_neighbor):(test_region_grid_coordinate(3)+num_neighbor));
% [idx_1, idx_2, idx_3] = ndgrid(3:5, 14:16, 58:60);
target_data_folder = DataManager.fp_raw_data_folder(dataset_name, stack);
for file_idx = 1 : numel(test_region_file_list)
    tic
    if isempty(test_region_file_list{file_idx})
        continue;
    end
    [test_folder_name, ~, ~]= fileparts(test_region_file_list{file_idx});
    target_folder_name = strrep(test_folder_name, '/nfs/birdstore-brainbucket2/Vessel/WholeBrain/mouselight_1/Raw_Green', target_data_folder);
    fprintf('Copying file %s\n', test_folder_name);
    mkdir(target_folder_name);
    copyfile(test_folder_name, target_folder_name);
    toc
end
disp('Finish copying files');
%% Load data
test_tile = DataManager.load_single_tiff(raw_data_grid.image_filepath{39});
test_tile = flip(test_tile, 1);
test_tile = flip(test_tile, 2);
test_tile2 = DataManager.load_single_tiff(raw_data_grid.image_filepath{40});
test_tile2 = flip(test_tile2, 1);
test_tile2 = flip(test_tile2, 2);
max_proj_1 = max(test_tile, [], 3);
max_proj_2 = max(test_tile2, [], 3);
figure;
imshow(max_proj_1)
figure;
imshow(max_proj_2)
est_disp = 79;
imshow(cat(2, max_proj_1(:, 1:(end - est_disp)), ...
    max_proj_2(:, (est_disp):end)))