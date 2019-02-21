%% Conclusion: Doesn't work. 
% This algorithm can be used to roughly compute the global translational
% transformation, but not for precious local translational transformation.
clc;clear;
DataManager = FileManager;
raw_data_info = load('~/Documents/Github/VesselReconstruction/Metadata/mouselight_1_raw_data_info.mat');
% dataset_folder = '/data/Vessel/ML_stitching/4_15_59_cube';
dataset_folder = '/data/Vessel/ML_stitching/4_17_57_cube';
raw_data_folder = fullfile(dataset_folder, 'raw_data');
descriptor_folder = fullfile(dataset_folder, 'stage_2_descriptor_output');
input_file_path = dir(fullfile(raw_data_folder, '**/*.tif'));
descriptor_filepaths = dir(fullfile(descriptor_folder, '**/*tor.mat'));
tile_idx_1 = 5;
tile_idx_2 = 14;
image_fp_1 = fullfile(input_file_path(tile_idx_1).folder, input_file_path(tile_idx_1).name);
image_fp_2 = fullfile(input_file_path(tile_idx_2).folder, input_file_path(tile_idx_2).name);
tile_descriptor_fp_1 = fullfile(descriptor_filepaths(tile_idx_1).folder, descriptor_filepaths(tile_idx_1).name);
tile_descriptor_fp_2 = fullfile(descriptor_filepaths(tile_idx_2).folder, descriptor_filepaths(tile_idx_2).name);
tile_acqusition_fd_1 = input_file_path(tile_idx_1).folder;
tile_acqusition_fd_2 = input_file_path(tile_idx_2).folder;
pixshift = [0,0,0];
image_1 = DataManager.load_single_tiff(image_fp_1);
image_2 = DataManager.load_single_tiff(image_fp_2);
desc_1 = load(tile_descriptor_fp_1);
desc_2 = load(tile_descriptor_fp_2);
%% Generate a 3D grid for the image
raw_image_size = [1536 1024 251];
image_size = raw_image_size([2,1,3]);
descriptor_1_sub = correctTiles(desc_1.blv_sub, image_size);
descriptor_2_sub = correctTiles(desc_2.blv_sub, image_size);
%% Load acquisition information to estimate the tile position
scopefile1 = readScopeFile(tile_acqusition_fd_1);
scopefile2 = readScopeFile(tile_acqusition_fd_2);
imsize_um = [scopefile1.x_size_um,scopefile1.y_size_um,scopefile1.z_size_um];
% estimate translation
gridshift = ([scopefile2.x scopefile2.y scopefile2.z]-[scopefile1.x scopefile1.y scopefile1.z]);
iadj =find(gridshift);
% Stage shift in micron
stgshift = 1000*([scopefile2.x_mm scopefile2.y_mm scopefile2.z_mm]-[scopefile1.x_mm scopefile1.y_mm scopefile1.z_mm]);
% Conver the stage shift to voxel distance
if all(pixshift==0)
    pixshift = round(stgshift.*(image_size-1)./imsize_um);
end


%% Transform the descriptor according to the estimated shift first
desc_2_sub_shifted = bsxfun(@plus, descriptor_2_sub, pixshift);
desc_2_sub_shifted_max = max(desc_2_sub_shifted, [], 1);
desc_2_sub_shifted_min = min(desc_2_sub_shifted, [], 1);

desc_1_sub_max = max(descriptor_1_sub, [], 1);
desc_1_sub_min = min(descriptor_1_sub, [], 1);
% Add offset to remove the boundary produced by cutting and region without
% signal
overlap_bbox_max = min(desc_2_sub_shifted_max, desc_1_sub_max) - 10;
overlap_bbox_min = max(desc_2_sub_shifted_min, desc_1_sub_min) + 10;
% overlap_bbox_min(1) = overlap_bbox_min(1) + 100;
desc_1_selected_Q = all(bsxfun(@ge, descriptor_1_sub, overlap_bbox_min) & bsxfun(@le, descriptor_1_sub, overlap_bbox_max), 2);
desc_2_selected_Q = all(bsxfun(@ge, desc_2_sub_shifted, overlap_bbox_min) & bsxfun(@le, desc_2_sub_shifted, overlap_bbox_max), 2);

desc_1_selected = descriptor_1_sub(desc_1_selected_Q,:);
desc_2_selected = desc_2_sub_shifted(desc_2_selected_Q,:);
%% Merge edge voxels
merge_box_size = [6,6,2];
desc_1_sub_ds = fun_stitching_merge_surface_voxels(desc_1_selected, merge_box_size);
desc_2_sub_ds = fun_stitching_merge_surface_voxels(desc_2_selected, merge_box_size);
% figure
% subplot(1,2,1)
% scatter3(desc_1_sub_ds(:,1), desc_1_sub_ds(:,2), desc_1_sub_ds(:,3))
% subplot(1,2,2)
% hold on
% scatter3(desc_2_sub_ds(:,1), desc_2_sub_ds(:,2), desc_2_sub_ds(:,3))
%% Test delany%
% Reconstruct the surfaces in 3D and triangulate it for each connected
% components
surface_mask_1 = false(image_size);
surface_mask_1(sub2ind(image_size, descriptor_1_sub(:,1), descriptor_1_sub(:,2), ...
    descriptor_1_sub(:,3))) = true;
surface_mask_1_cc = bwconncomp(surface_mask_1);
test_cc_sub = fun_ind2sub(image_size, surface_mask_1_cc.PixelIdxList{3});
surface_mask_1_cc_tri = delaunayTriangulation(test_cc_sub(:,1), test_cc_sub(:,2), test_cc_sub(:,3));

%% Non-regid matching with low rank approximation
projectionThr = 5;
optimopts = statset('nlinfit');
optimopts.RobustWgtFun = 'bisquare';
opt.method='nonrigid_lowrank';
opt.beta=6;            % the width of Gaussian kernel (smoothness)
opt.lambda=16;          % regularization weight
opt.viz=0;              % show every iteration
opt.outliers=0.9;       % use 0.7 noise weight
opt.numeig = 100;   
opt.eigfgt = 0;
opt.max_it = 50;
opt.fgt=0;              % do not use FGT (default)
opt.tol = 1e-3;
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)
matchparams.optimopts = optimopts;
matchparams.opt = opt;
matchparams.projectionThr = projectionThr;
matchparams.model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
matchparams.debug = false;
matchparams.viz = false;
tic
[rate,X_,Y_,tY_] = vessel_descriptorMatchforz(X,Y,pixshift,iadj,matchparams);
toc
% save('~/Documents/Github/MouseLight/Test_stitching/surface_registration_test_nonrigid.mat', 'desc_1_sub_ds', 'desc_2_sub_ds', ...
%     'X', 'Y', 'X_', 'Y_', 'tY_', 'rate');
tmp = load('~/Documents/Github/MouseLight/Test_stitching/surface_registration_test_nonrigid.mat', 'desc_1_sub_ds', 'desc_2_sub_ds', ...
    'X', 'Y', 'X_', 'Y_', 'tY_', 'rate');

% scatter3(desc_1_sub_ds(:,1), desc_1_sub_ds(:,2), desc_1_sub_ds(:,3))
% hold on 
% scatter3(desc_2_sub_ds(:,1), desc_2_sub_ds(:,2), desc_2_sub_ds(:,3))
scatter3(X(:,1), X(:,2), X(:,3))
hold on 
scatter3(Y(:,1), Y(:,2), Y(:,3))
xlabel 'x'
ylabel 'y'
zlabel 'z'
legend('Tile 1', 'Tile 2');

figure
scatter3(X_lr(:,1), X_lr(:,2), X_lr(:,3))
hold on 
scatter3(tY_lr(:,1), tY_lr(:,2), tY_lr(:,3))
xlabel 'x'
ylabel 'y'
zlabel 'z'
legend('Tile 1', 'Tile 2');

% Displacement
figure;
subplot(1,3,1);
histogram(X_lr(:,1) - Y_lr(:,1))
title('Displacement in X')
subplot(1,3,2);
histogram(X_lr(:,2) - Y_lr(:,2))
title('Displacement in Y')
subplot(1,3,3);
histogram(X_lr(:,3) - Y_lr(:,3))
title('Displacement in Z')
num_edge_voxel = zeros(grid_3d.size);
mean_edge_pos_1 = zeros(grid_3d.size);
mean_edge_pos_2 = zeros(grid_3d.size);
mean_edge_pos_3 = zeros(grid_3d.size);
disp_field = X_ - Y_;
bbox_sub = ceil(X_ ./ block_size);
for iter1 = 1 : size(bbox_sub, 1)
    num_edge_voxel(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
        num_edge_voxel(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + 1;
    mean_edge_pos_1(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
        mean_edge_pos_1(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + disp_field(iter1, 1);
    mean_edge_pos_2(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
        mean_edge_pos_2(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + disp_field(iter1, 2);
    mean_edge_pos_3(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
        mean_edge_pos_3(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + disp_field(iter1, 3);
end
vis_sec = 50;
figure;
subplot(1,3,1)
imagesc(mean_edge_pos_1(:,:,vis_sec))
subplot(1,3,2)
imagesc(mean_edge_pos_2(:,:,vis_sec))
subplot(1,3,3)
imagesc(mean_edge_pos_3(:,:,vis_sec))
% tmp = load('~/Documents/Github/MouseLight/Test_stitching/surface_registration_test_nonrigid.mat');
disp_field = X_lr - Y_lr;
median_disp = mode(disp_field,1);
quiver3(X_lr(:,1), X_lr(:,2), X_lr(:,3), disp_field(:,1) - median_disp(1), ...
    disp_field(:,2) - median_disp(2), disp_field(:,3) - median_disp(3), 'AutoScaleFactor', 2);