clc;clear;
DataManager = FileManager;
raw_data_info = load('~/Documents/Github/VesselReconstruction/Metadata/mouselight_1_raw_data_info.mat');
dataset_folder = '/data/Vessel/ML_stitching/4_15_59_cube';
raw_data_folder = '/data/Vessel/ML_stitching/4_15_59_cube/raw_data';
descriptor_folder = fullfile(dataset_folder, 'stage_2_descriptor_output_ori');
input_file_path = dir(fullfile(raw_data_folder, '**/*.tif'));
scope_files = dir(fullfile(raw_data_folder, '**/*.acquisition'));
scope_file_names = cellfun(@fullfile, {scope_files.folder}, {scope_files.name}, 'UniformOutput', false);
% Output scope acquisition list ( required by pointmatch_task
% fid = fopen(fullfile(raw_data_folder, 'scopeacquisitionlist.txt'), 'w');
% for idx = 1 : numel(scope_file_names)
%     fprintf(fid, '%s\n', scope_file_names{idx});
% end
% fclose(fid);

num_file = numel(input_file_path);
configfile = '/home/dklab/Documents/Github/MouseLight/pipeline-descriptor-master/configfiles/2018-08-15.cfg';
%% Generate descriptor

for file_idx = [5, 14, 23]
    output_folder = strrep(input_file_path(file_idx).folder, raw_data_folder, descriptor_folder);
    output_filename = strrep(input_file_path(file_idx).name, '.tif', 'descriptor.mat');
    image_fp = fullfile(input_file_path(file_idx).folder, input_file_path(file_idx).name);
    fprintf('Processing %s\n', image_fp);
    output_fp = fullfile(output_folder, output_filename);
    if ~isfolder(output_folder)
        mkdir(output_folder);
    end
    tic
%     exit_code = vesselDescriptor_v2(image_fp, output_fp, configfile);
    exit_code = skelDescriptor(image_fp, output_fp, configfile);
    toc
end

%% Point Match - between two tiles
matching_folder = fullfile(dataset_folder, 'stage_3_point_match_output');
descriptor_filepaths = dir(fullfile(descriptor_folder, '**/*tor.txt'));
tile_idx_1 = 14;
tile_idx_2 = 23;
tile_descriptor_fp_1 = fullfile(descriptor_filepaths(2).folder, descriptor_filepaths(2).name);
tile_descriptor_fp_2 = fullfile(descriptor_filepaths(3).folder, descriptor_filepaths(3).name);
tile_image_fp_1 = fullfile(input_file_path(tile_idx_1).folder, input_file_path(tile_idx_1).name);
tile_image_fp_2 = fullfile(input_file_path(tile_idx_2).folder, input_file_path(tile_idx_2).name);
%% Coordinate check 
% tile_image_1 = DataManager.load_single_tiff(tile_image_fp_1);
% tile_image_2 = DataManager.load_single_tiff(tile_image_fp_2);
% tile_image_maxpro_flip_1 = rot90(max(tile_image_1, [], 3),2);
% tile_image_maxpro_flip_2 = rot90(max(tile_image_2, [], 3),2);
% tile_size = size(tile_image_1);
% desc1 = readDesc(tile_descriptor_fp_1,{'0'});
% desc2 = readDesc(tile_descriptor_fp_2,{'0'});
% desc1 = correctTiles(desc1,tile_size);
% desc2 = correctTiles(desc2,tile_size);
% 
% subplot(1,4,1)
% imshow(tile_image_maxpro_flip_1);
% subplot(1,4,2);
% imshow(tile_image_maxpro_flip_2);
% subplot(1,4,3);
% scatter(desc1(:,2), desc1(:,1));
% set(gca, 'YDir', 'reverse');
% subplot(1,4,4);
% scatter(desc2(:,2), desc2(:,1));
% set(gca, 'YDir', 'reverse');
% ind_2 = sub2ind(tile_size, desc2(:,1)+1, desc2(:,2)+1, desc2(:,3)+1);
% skel_mask = false(tile_size);
% skel_mask(ind_2) = true;
% DataManager.visualize_itksnap(tile_image_2, skel_mask);
%%

tile_acqusition_fd_1 = input_file_path(tile_idx_1).folder;
tile_acqusition_fd_2 = input_file_path(tile_idx_2).folder;
output_image_folder = strrep(input_file_path(tile_idx_1).folder, raw_data_folder, matching_folder);
if ~isfolder(output_image_folder)
    mkdir(output_image_folder);
end
descriptor_channel = {'0'};
pixshift = [0,0,0];
% profile on 
[~] = pointmatch(tile_descriptor_fp_1, tile_descriptor_fp_2, tile_acqusition_fd_1, tile_acqusition_fd_2, ...
    output_image_folder,pixshift, descriptor_channel);
% [~] = pointmatch_vessel(tile_descriptor_fp_1, tile_descriptor_fp_2, tile_acqusition_fd_1, tile_acqusition_fd_2, ...
%     output_image_folder,pixshift, descriptor_channel);
% profile off
% profile viewer
%% Visualization matching 
scatter3(X_(:,1), X_(:,2), X_(:,3))
hold on
scatter3(Y_(:,1) + pixshiftout(1), Y_(:,2) + pixshiftout(2), Y_(:,3) + pixshiftout(3))