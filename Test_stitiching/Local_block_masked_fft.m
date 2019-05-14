%% Conclusion: Doesn't work. 
% This algorithm can be used to roughly compute the global translational
% transformation, but not for precious local translational transformation.
DataManager = FileManager;
tile_idx_1 = 5;
tile_idx_2 = 14;
image_fp_1 = fullfile(input_file_path(tile_idx_1).folder, input_file_path(tile_idx_1).name);
image_fp_2 = fullfile(input_file_path(tile_idx_2).folder, input_file_path(tile_idx_2).name);
tile_descriptor_fp_1 = fullfile(descriptor_filepaths(tile_idx_1).folder, descriptor_filepaths(tile_idx_1).name);
tile_descriptor_fp_2 = fullfile(descriptor_filepaths(tile_idx_2).folder, descriptor_filepaths(tile_idx_2).name);
image_1 = DataManager.load_single_tiff(image_fp_1);
image_2 = DataManager.load_single_tiff(image_fp_2);
desc_1 = load(tile_descriptor_fp_1);
desc_2 = load(tile_descriptor_fp_2);
%% Generate a 3D grid for the image
image_size = size(image_1);
block_size = 32;
search_pad_size = 16;
grid_3d = fun_generate_grid(block_size, 0, image_size);
%% Determine the position of the boundary
% 1. Use edge output from previous descriptor
% 2. Later, find the boundary of the large vessel mask
tile_displacement = [0, 0, 155];
grid_3d.num_edge_voxel = zeros(grid_3d.size);
% Select the descriptor position in tile 1 in the overlapping region
edge_sub_kept_Q = all(bsxfun(@gt, desc_1.blv_sub, tile_displacement), 2);
edge_sub_kept = desc_1.blv_sub(edge_sub_kept_Q,[2,1,3]);
mask_fft_block_size = [block_size,block_size,block_size];
% Divide the data in the overlapping region from tile 1 into boxes

bbox_sub_1 = ceil(edge_sub_kept(:,1) / mask_fft_block_size(1));
bbox_sub_2 = ceil(edge_sub_kept(:,2) / mask_fft_block_size(2));
bbox_sub_3 = ceil(edge_sub_kept(:,3) / mask_fft_block_size(3));
% Select block for matching
for iter1 = 1 : numel(bbox_sub_1)
    grid_3d.num_edge_voxel(bbox_sub_1(iter1), bbox_sub_2(iter1), bbox_sub_3(iter1)) = ...
        grid_3d.num_edge_voxel(bbox_sub_1(iter1), bbox_sub_2(iter1), bbox_sub_3(iter1)) + 1;
end
grid_3d.valid_grid_ind = find(grid_3d.num_edge_voxel > 0);
grid_3d.valid_grid_num_edge_voxel = grid_3d.num_edge_voxel(grid_3d.valid_grid_ind);
[grid_3d.valid_grid_num_edge_voxel, tmp_idx] = sort(grid_3d.valid_grid_num_edge_voxel, 'descend');
grid_3d.valid_grid_ind = grid_3d.valid_grid_ind(tmp_idx(grid_3d.valid_grid_num_edge_voxel > 1100));
num_valid_grid = numel(grid_3d.valid_grid_ind);
scan_shift = zeros(3, num_valid_grid);
scan_score = zeros(1, num_valid_grid);
tic
for test_grid_idx = 1 : num_valid_grid 
    test_grid_ind = grid_3d.valid_grid_ind(test_grid_idx);
    test_grid_pos = grid_3d.center_pos(test_grid_ind, :);
    test_bbox_1 = grid_3d.mmll(test_grid_ind,:);
    test_block_image = crop_bbox3(image_1, test_bbox_1, 'default');
    bbox_tile_2_mmxx = grid_3d.mmxx(test_grid_ind,:);
    bbox_tile_2_mmxx(1:3) = max(1, bbox_tile_2_mmxx(1:3) - search_pad_size - tile_displacement);
    bbox_tile_2_mmxx(4:6) = min(image_size, bbox_tile_2_mmxx(4:6) + search_pad_size - tile_displacement);
    bbox_tile_2_mmll = bbox_tile_2_mmxx;
    bbox_tile_2_mmll(4:6) = bbox_tile_2_mmll(4:6) - bbox_tile_2_mmll(1:3);
    bbox_displacement = bbox_tile_2_mmll(1:3) - test_bbox_1(1:3) + tile_displacement;
    test_target_image = crop_bbox3(image_2, bbox_tile_2_mmll, 'default');
    fixed_mask = test_target_image > 1.5e4;
    moving_mask = test_block_image > 1.5e4;
    
    [transform, scan_score(test_grid_idx), ~, ~] = MaskedTranslationRegistration(test_target_image,test_block_image,fixed_mask,moving_mask);
    scan_shift(:, test_grid_idx) = (transform + bbox_displacement)';
end
toc
%% Select high score shift
high_score_Q = scan_score > 0.9;
shift_x = scan_shift(1, high_score_Q);
shift_y = scan_shift(2, high_score_Q);
shift_z = scan_shift(3, high_score_Q);
%% Visualize displace field
score_field = zeros(grid_3d.size);
dis_field_x = ones(grid_3d.size) .* 0;
dis_field_y = ones(grid_3d.size) .* 0;
dis_field_z = ones(grid_3d.size) .* 0;
score_field(grid_3d.valid_grid_ind) = scan_score;
dis_field_x(grid_3d.valid_grid_ind) = scan_shift(1,:) - mode(scan_shift(1,high_score_Q));
dis_field_y(grid_3d.valid_grid_ind) = scan_shift(2,:) - mode(scan_shift(2,high_score_Q));
dis_field_z(grid_3d.valid_grid_ind) = scan_shift(3,:) - mode(scan_shift(3,high_score_Q));
%%
layer = 8;
score_th = 0.9;
figure;
subplot(1,4,1);
imagesc(score_field(:,:,layer).* (score_field(:,:,layer) > score_th));
colormap jet
colorbar
subplot(1,4,2);
imagesc(dis_field_x(:,:,layer).* (score_field(:,:,layer) > score_th));
colorbar
subplot(1,4,3);
imagesc(dis_field_y(:,:,layer).* (score_field(:,:,layer) > score_th));
colorbar
subplot(1,4,4);
imagesc(dis_field_z(:,:,layer).* (score_field(:,:,layer) > score_th));
colorbar


% overlapping_image_size = max(0, image_size - tile_displacement);


% scatter3(edge_sub_kept(:,1), edge_sub_kept(:,2), edge_sub_kept(:,3))

