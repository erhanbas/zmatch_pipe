raw_data_root = '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2018-08-15';
% Gather raw images file paths
raw_data_info = dir(fullfile(raw_data_root, '**/*.0.tif'));
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
is_duplicate_file_Q = false(1, numel(raw_data_info));
for tile_idx = 1 : numel(raw_data_info)
    if contains(raw_data_info(tile_idx).folder, 'Duplicates')
        is_duplicate_file_Q(tile_idx) = true;
    end
end
raw_data_info = raw_data_info(~is_duplicate_file_Q);   
num_file = numel(raw_data_info);
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
save('/groups/mousebrainmicro/home/jix/Documents/GitHub/pipeline-featmatch/Test_stitiching/mouselight_1_raw_data_info.mat', '-struct', 'raw_data_grid');
%% Determine region of interest for testing stitching algorithm
% Region with both large vessels and capillaries: 
clc;clear;
raw_data_grid = load('/groups/mousebrainmicro/home/jix/Documents/GitHub/pipeline-featmatch/Test_stitiching/mouselight_1_raw_data_info.mat');
% DataManager = FileManager;
% dataset_name = 'ML_stitching';
% for iter_tile = 1 : raw_data_grid.num_tile
%     if ~isempty(strfind(raw_data_grid.image_filepath{iter_tile}, '2018-08-23/01/01567'))
%         tile_idx = iter_tile;
%         test_region_grid_coordinate = raw_data_grid.grid_pos_sub(iter_tile, :);
%     end
% end

% for iter_tile = 1 : prod(raw_data_grid.grid_size)
%     if ~isempty(strfind(raw_data_grid.image_filepath_array{iter_tile}, '2018-08-23/01/01567'))
%         tile_ind = iter_tile;
% %         test_region_grid_coordinate = raw_data_grid.grid_pos_sub(:, iter_tile)';
%     end
% end
test_region_grid_coordinate = [6, 19, 58];
stack = [num2str(test_region_grid_coordinate, '%02d_'), 'cube'];
num_neighbor = 3;
test_region_file_list = raw_data_grid.image_filepath_array((test_region_grid_coordinate(1)-num_neighbor):(test_region_grid_coordinate(1)+num_neighbor),...
    (test_region_grid_coordinate(2)-num_neighbor):(test_region_grid_coordinate(2)+num_neighbor), ...
    (test_region_grid_coordinate(3)-num_neighbor):(test_region_grid_coordinate(3)+num_neighbor));
% [idx_1, idx_2, idx_3] = ndgrid(3:5, 14:16, 58:60);
% target_data_folder = DataManager.fp_raw_data_folder(dataset_name, stack);
target_data_folder = fullfile('/nrs/mouselight/Users/jix/pipeline_test', stack, 'raw_data');
for file_idx = 1 : numel(test_region_file_list)
    tic
    if isempty(test_region_file_list{file_idx})
        continue;
    end
    source_tiff_fp = test_region_file_list{file_idx};
    source_acquisition_fp = strrep(source_tiff_fp, '0.tif', 'acquisition');
    source_microscope_fp = strrep(source_tiff_fp, '0.tif', 'microscope');
    target_tiff_fp = strrep(source_tiff_fp, '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2018-08-15', target_data_folder);
    target_acquisition_fp = strrep(source_acquisition_fp, '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2018-08-15', target_data_folder);
    target_microscope_fp = strrep(source_microscope_fp, '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2018-08-15', target_data_folder);
    target_folder_name = fileparts(target_tiff_fp);

    fprintf('Copying file %s\n', test_folder_name);
    mkdir(target_folder_name);
    copyfile(source_tiff_fp, target_tiff_fp);
    copyfile(source_acquisition_fp, target_acquisition_fp);
    copyfile(source_microscope_fp, target_microscope_fp);
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