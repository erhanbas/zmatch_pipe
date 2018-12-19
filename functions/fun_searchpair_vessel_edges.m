function [X_stable,Y_stable,rate, pixshift] = fun_searchpair_vessel_edges(descriptor_1_sub,descriptor_2_sub,pixshiftinit)
% fun_searchpair_vessel_edges apply Coherent Point Drift for two point
% cloud (surface) registration. Since surface point cloud is dense, the
% point clouds are downsampled by mergering nearby points and reolaced by
% their mean position( fun_stitching_merge)surface_voxels ). The resulting
% point clouds are fed into CPD. The output voxel pairs are further
% selected by compute the local displacement standard deviation. Only local
% blocks with standard deviation smaller than a threshold are kept. 
% Input: 
%   descriptor_1_sub: N-by-3 numerical array, position of the points
%   descriptor_2_sub: M-by-3 numerical array, position of the points in the
%   other point cloud. 
%   pixshiftinit: 3-by-1 numerical array, initial estimation of the
%   translation displacment between two point clouds
% Output: 
%   X_stable: N-by-3 numerical array, position of the matched points
%   Y_stable: N-by-3 numerical array, position of the matched points
%   rate: numerical scalar, defined by Erhan, hard to explain... see
%   vessel_descriptorMatchforz. A value that reveals the reliable of
%   matching ( the higher, the better ) 
%   pixshift: 3-by-1 numerical array, median of the displacement of the
%   output paired point sets. 
% 
% Author: Xiang Ji ( xiangji.ucsd@gmail.com )
% Date: Dec 5, 2018

% If the number of available feature points is larger than this number,
% downsample the point cloud.
th_num_ds_desc = 30000; 

pixshift =  pixshiftinit;
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
desc_1_selected_Q = all(bsxfun(@ge, descriptor_1_sub, overlap_bbox_min) & bsxfun(@le, descriptor_1_sub, overlap_bbox_max), 2);
desc_2_selected_Q = all(bsxfun(@ge, desc_2_sub_shifted, overlap_bbox_min) & bsxfun(@le, desc_2_sub_shifted, overlap_bbox_max), 2);

desc_1_selected = descriptor_1_sub(desc_1_selected_Q,:);
desc_2_selected = desc_2_sub_shifted(desc_2_selected_Q,:);
%% Merge edge voxels
merge_box_size = [6, 6, 2];
if (size(desc_1_selected,1) > th_num_ds_desc ) || (size(desc_2_selected,1) > th_num_ds_desc)
    desc_1_sub_ds = fun_stitching_merge_surface_voxels(desc_1_selected, merge_box_size);
    desc_2_sub_ds = fun_stitching_merge_surface_voxels(desc_2_selected, merge_box_size);
else
    desc_1_sub_ds = desc_1_selected;
    desc_2_sub_ds = desc_2_selected;
end
if isempty(desc_1_sub_ds) || isempty(desc_2_sub_ds)
    X_stable = [];
    Y_stable = [];
    rate = 0;
    return;
end    
%% Non-regid matching with low rank approximation
projectionThr = 5;
optimopts = statset('nlinfit');
optimopts.RobustWgtFun = 'bisquare';
opt.method='nonrigid_lowrank';
opt.beta=6;            % the width of Gaussian kernel (smoothness)
opt.lambda=16;          % regularization weight
opt.viz=0;              % show every iteration
opt.outliers=0.9;       % use 0.7 noise weight
opt.numeig = 100;       % Number of eigenvectors for low rank appriximation 
opt.eigfgt = 0;         % Use fast gaussian transformation for computing eigenvectors
opt.max_it = 50;        % Maxinum number of iteration 
opt.fgt=0;              % do not use FGT (default)
opt.tol = 1e-3;         % Error tolorance ( how is it defined? )
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)
matchparams.optimopts = optimopts;
matchparams.opt = opt;
matchparams.projectionThr = projectionThr;
matchparams.model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
matchparams.debug = false;
matchparams.viz = false;
tic
[rate, X_, Y_, tY_] = vessel_descriptorMatchforz(desc_1_sub_ds, desc_2_sub_ds, pixshift, matchparams);
toc
%% Matched point selection 
%  If the matching is bad, return empay matching directly
if rate < 0.8 || isempty(X_) || isempty(Y_)
    X_stable = [];
    Y_stable = [];
    rate = 0;
    return;
end
% Check the consistancy of the local transformation 
% Compute the standard deviation of the displacement of the matched pairs
% in small blocks ( say, of size 10x10x10 um ) 
% displacement between X_ and Y_ for each local block 
disp_X_Y = X_ - Y_;
sub_min = min(X_, [], 1);
sub_max = max(X_, [], 1);
sub_bbox = sub_max - sub_min + 1;
local_bbox_size = [30,30,10];
grid_size = ceil(sub_bbox ./ local_bbox_size);
grid_sub = ceil((1 + bsxfun(@minus, X_, sub_min))./ local_bbox_size);
grid_ind = sub2ind(grid_size, grid_sub(:,1), grid_sub(:,2), grid_sub(:,3));
[grid_idx_cell, grid_idx_list]= fun_bin_data_to_idx_list(grid_ind);
num_valid_bbox = numel(grid_idx_list);
grid_idx_cell_array = cell(grid_size);
num_pair = zeros(grid_size);
% mean_dis_1 = zeros(grid_size);
% mean_dis_2 = zeros(grid_size);
% mean_dis_3 = zeros(grid_size);
std_dis_1 = zeros(grid_size) - 3;
std_dis_2 = zeros(grid_size) - 3;
std_dis_3 = zeros(grid_size) - 3;
for iter_bbox = 1 : num_valid_bbox
    dis_list = disp_X_Y(grid_idx_cell{iter_bbox},:);
    tmp_grid_ind = grid_idx_list(iter_bbox);
    tmp_num_pair = size(dis_list, 1);
    num_pair(tmp_grid_ind) = tmp_num_pair;
    grid_idx_cell_array{tmp_grid_ind} = grid_idx_cell{iter_bbox};
    if tmp_num_pair > 1
%         mean_dis_1(tmp_grid_ind) = mean(dis_list(:,1));
%         mean_dis_2(tmp_grid_ind) = mean(dis_list(:,2));
%         mean_dis_3(tmp_grid_ind) = mean(dis_list(:,3));
        std_dis_1(tmp_grid_ind) = std(dis_list(:,1));
        std_dis_2(tmp_grid_ind) = std(dis_list(:,2));
        std_dis_3(tmp_grid_ind) = std(dis_list(:,3));
    elseif tmp_num_pair == 1
%         mean_dis_1(tmp_grid_ind) = dis_list(1);
%         mean_dis_2(tmp_grid_ind) = dis_list(:,2);
%         mean_dis_3(tmp_grid_ind) = dis_list(:,3);
    end
end
% Select the matches in the low standard deviation grids ( roughly less
% than 0.5 micron) 
grid_low_std_Q = (std_dis_1 <= 1.5) & (std_dis_2 <= 1.5) & (std_dis_3 <= 0.5) & ...
    ( num_pair > 1 );
grid_low_std_voxel_idx = cat(2, grid_idx_cell_array{grid_low_std_Q});
X_stable = X_(grid_low_std_voxel_idx,:);
Y_stable = Y_(grid_low_std_voxel_idx,:);
% The initial estimation of the stage position is quite close. If the
% displacement between Y_ and tY_ is too large, it must be a wrong
% matching
disp_X_Y = X_stable - Y_stable;
disp_X_Y_med = median(disp_X_Y);
disp_X_Y_dev = disp_X_Y - disp_X_Y_med;
disp_X_Y_dev_std = std(single(disp_X_Y_dev),1);
disp_inlier = all(abs(disp_X_Y_dev) < disp_X_Y_dev_std .* 2, 2);
X_stable = X_stable(disp_inlier);
Y_stable = Y_stable(disp_inlier);
pixshift = median(X_stable - Y_stable);
%% Visualize matched points
% figure;
% subplot(1,3,1)
% scatter3(desc_1_selected(:, 1), desc_1_selected(:, 2), desc_1_selected(:, 3));
% hold on
% scatter3(desc_2_selected(:, 1), desc_2_selected(:, 2), desc_2_selected(:, 3));
% xlabel 'x';
% ylabel 'y';
% zlabel 'z';
% legend('Tile 1', 'Tile 2');
% title('Descriptors in the overlapping region');
% subplot(1,3,2)
% scatter3(X_(:, 1), X_(:, 2), X_(:, 3));
% hold on
% scatter3(Y_(:, 1) + pixshiftinit(1), Y_(:, 2) + pixshiftinit(2), Y_(:, 3) + pixshiftinit(3));
% xlabel 'x';
% ylabel 'y';
% zlabel 'z';
% legend('Tile 1', 'Tile 2');
% title('Matched Descriptors in the overlapping region');
% subplot(1,3,3)
% scatter3(X_stable(:,1), X_stable(:,2), X_stable(:,3))
% hold on 
% scatter3(Y_stable(:,1) + pixshiftinit(1), Y_stable(:,2) + pixshiftinit(2), Y_stable(:,3) + pixshiftinit(3))
% legend('Tile 1', 'Tile 2');
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('Matched edge voxel with low local variance in displacement')
%% Visualization of the displacement field
% vis_sec = 1;
% figure;
% subplot(2,3,1)
% imagesc(mean_dis_1(:,:,vis_sec));
% title('Mean displacement in X')
% colorbar
% subplot(2,3,2)
% imagesc(mean_dis_2(:,:,vis_sec));
% title('Mean displacement in Y')
% colorbar
% subplot(2,3,3)
% imagesc(mean_dis_3(:,:,vis_sec));
% title('Mean displacement in Z')
% colorbar
% subplot(2,3,4)
% imagesc(std_dis_1(:,:,vis_sec));
% title('STD of displacement in X');
% colorbar
% subplot(2,3,5)
% imagesc(std_dis_2(:,:,vis_sec));
% title('STD of displacement in Y')
% colorbar
% subplot(2,3,6)
% imagesc(std_dis_3(:,:,vis_sec));
% title('STD of displacement in Z')
% colorbar
% figure;
% scatter3(X_stable(:,1), X_stable(:,2), X_stable(:,3))
% hold on 
% scatter3(Y_stable(:,1), Y_stable(:,2), Y_stable(:,3))
% legend('Tile 1', 'Tile 2');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('Matched edge voxel with low local variance in displacement')
% view(2)
%% 
end

%% Sub function
function voxel_sub = fun_stitching_merge_surface_voxels(voxel_sub, merge_block_size)
% fun_stitching_merge_surface_voxels merges the input voxel subscript list
% by computing the average position of the voxels within the block, whose
% size is specified by merge_block_size
% Input: 
%   voxel_sub: N-by-3, coordinate of the voxel position in 3D space
%   merge_block_size: 1-by-3 numerical vector, size of the block for merging. 
% Output: 
%   voxel_sub: N'-by-3 numerical array after merging
% Author: Xiang Ji (xiangji.ucsd@gmail.com) 
% Date: Dec 5, 2018
sub_min = min(voxel_sub, [], 1);
sub_max = max(voxel_sub, [], 1);
image_size = sub_max - sub_min + 1;
downsampled_image_size = ceil(image_size ./ merge_block_size);

num_edge_voxel = zeros(downsampled_image_size);
mean_edge_pos_1 = zeros(downsampled_image_size);
mean_edge_pos_2 = zeros(downsampled_image_size);
mean_edge_pos_3 = zeros(downsampled_image_size);
bbox_sub = ceil((1 + bsxfun(@minus, voxel_sub, sub_min))./ merge_block_size);
for iter1 = 1 : size(bbox_sub, 1)
    num_edge_voxel(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
        num_edge_voxel(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + 1;
    mean_edge_pos_1(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
        mean_edge_pos_1(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + voxel_sub(iter1, 1);
    mean_edge_pos_2(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
        mean_edge_pos_2(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + voxel_sub(iter1, 2);
    mean_edge_pos_3(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
        mean_edge_pos_3(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + voxel_sub(iter1, 3);
end
voxel_sub = cat(2, mean_edge_pos_1(:) ./ max(1, num_edge_voxel(:)), ...
    mean_edge_pos_2(:) ./ max(1, num_edge_voxel(:)), mean_edge_pos_3(:) ./ max(1, num_edge_voxel(:)));

voxel_sub = voxel_sub(all(voxel_sub > 0, 2),:);
end
%% Sub function
function [rate,X_,Y_,tY_] = vessel_descriptorMatchforz(X,Y,pixshift,params)
%DESCRIPTORMATCH Summary of this function goes here
%
% [OUTPUTARGS] = DESCRIPTORMATCH(INPUTARGS) Explain usage here
%
% Inputs:
%   X, Y: two 2D real, double marices, specifying the position of the point
%   in the point cloud.
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2016/09/23 14:09:29 $	$Revision: 0.1 $
% Copyright: HHMI 2016
opt = params.opt;
projectionThr = params.projectionThr;
%% Initial match based on point drift
[Transform, ~] = cpd_register(X,Y,opt);
%% check if match is found
% Compute the pairwise euclidean distance between two input array
pD = pdist2(X,Transform.Y);
[aa1,bb1] = min(pD,[],1);
[~,bb2] = min(pD,[],2);
keeptheseY = find([1:length(bb1)]'==bb2(bb1));
keeptheseX = bb1(keeptheseY)';

disttrim = aa1(keeptheseY)' < projectionThr;
X_ = X(keeptheseX(disttrim),:);
Y_ = Y(keeptheseY(disttrim),:);
tY_= Transform.Y(keeptheseY(disttrim),:);
% Rate is the ratio of the number of matched pair of distance less than
% projectionThr over the total number of matched pairs
rate = sum(disttrim)/length(disttrim);
% [pixshift rate]
if rate < .5 % dont need to continue
    [X_,Y_,~] = deal(0);
    return
end
Y_ = bsxfun(@minus, Y_, pixshift);
end