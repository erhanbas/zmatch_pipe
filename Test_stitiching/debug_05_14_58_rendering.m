clc;clear;
Vessel_data_root_path = '/data/Vessel/';
% raw_data_info = load(fullfile(DataManager.SCRIPT_PATH, 'Metadata',  'mouselight_1_raw_data_info.mat'));
brain = '05_14_58_cube';
dataset_folder = fullfile(Vessel_data_root_path, 'ML_stitching', brain);
raw_data_folder = fullfile(dataset_folder, 'raw_data');
descriptor_folder = fullfile(dataset_folder, 'stage_2_descriptor_output');
input_file_path = dir(fullfile(raw_data_folder, ['**', filesep, '*.tif']));
scope_files = dir(fullfile(raw_data_folder, ['**', filesep, '*.acquisition']));
scope_file_names = cellfun(@fullfile, {scope_files.folder}, {scope_files.name}, 'UniformOutput', false);
% Output scope acquisition list ( required by pointmatch_task
% fid = fopen(fullfile(raw_data_folder, 'scopeacquisitionlist.txt'), 'w');
% for idx = 1 : numel(scope_file_names)
%     fprintf(fid, '%s\n', scope_file_names{idx});
% end
% fclose(fid);

num_file = numel(input_file_path);
configfile = '/home/dklab/Documents/Github/MouseLight/pipeline-descriptor-master/configfiles/2018-08-15.cfg';
% Load raw data information
mouselight_raw_data_info = load('~/Documents/Github/VesselReconstruction/Metadata/mouselight_1_raw_data_info.mat');
mouselight_data_info = load('~/Documents/Github/VesselReconstruction/Metadata/mouselight_1_data_info.mat');
mouselight_file_str = load('~/Documents/Github/VesselReconstruction/Metadata/mouselight_file_str.mat');
%% 
check_tile = [5,7,7,2];
chekc_tile_array_sub = fun_io_octree_coordinate_to_array_sub(check_tile);
tile_size = [576 760 72];
check_tile_pix_sub_orgin = (chekc_tile_array_sub - 1) .* tile_size;
check_tile_pix_sub = check_tile_pix_sub_orgin + [280 420 34];
render_voxel_size = [0.24887814670138889,0.24949786184210523,1.0323237847222222];
raw_voxel_size = mouselight_raw_data_info.voxel_size;

render_data_size = tile_size .* 16;
render_data_size_um = render_data_size .* render_voxel_size;
stitched_raw_data_size = round(render_data_size_um ./ raw_voxel_size);

% est_raw_tile_array_sub = (check_tile_pix_sub ./ stitched_raw_data_size) * 7;
est_raw_tile_array_sub = [7,7,7];
est_raw_tile_global_array_sub = [5 14 58] - 4 + est_raw_tile_array_sub;
mouselight_raw_data_info.image_filepath_array{est_raw_tile_global_array_sub(1), ...
    est_raw_tile_global_array_sub(2), est_raw_tile_global_array_sub(3)}


mouselight_raw_data_info.image_filepath_array{5,15,57}
%% Generate descriptor
for file_idx = 1 : num_file
    output_folder = strrep(input_file_path(file_idx).folder, raw_data_folder, descriptor_folder);
    output_filename = strrep(input_file_path(file_idx).name, '.tif', 'descriptor.mat');
    image_fp = fullfile(input_file_path(file_idx).folder, input_file_path(file_idx).name);
    fprintf('Processing %s\n', image_fp);
    output_fp = fullfile(output_folder, output_filename);
    if ~isfolder(output_folder)
        mkdir(output_folder);
    end
    tic
    exit_code = vesselDescriptor(image_fp, output_fp, configfile);
%     exit_code = skelDescriptor(image_fp, output_fp, configfile);
    toc
end
%% Point Match - between two tiles
matching_folder = fullfile(dataset_folder, 'stage_3_point_match_output');
descriptor_filepaths = dir(fullfile(descriptor_folder, '**', '*tor.mat'));
tile_idx_1 = 3;
tile_idx_2 = 12;
tile_descriptor_fp_1 = fullfile(descriptor_filepaths(tile_idx_1).folder, descriptor_filepaths(tile_idx_1).name);
tile_descriptor_fp_2 = fullfile(descriptor_filepaths(tile_idx_2).folder, descriptor_filepaths(tile_idx_2).name);
tile_acqusition_fd_1 = scope_files(tile_idx_1).folder;
tile_acqusition_fd_2 = scope_files(tile_idx_2).folder;
output_image_folder = strrep(scope_files(tile_idx_1).folder, raw_data_folder, matching_folder);
if ~isfolder(output_image_folder)
    mkdir(output_image_folder);
end
descriptor_channel = {'0'};
pixshift = [0,0,0];
[~] = pointmatch_vessel(tile_descriptor_fp_1, tile_descriptor_fp_2, tile_acqusition_fd_1, tile_acqusition_fd_2, ...
    output_image_folder,pixshift, descriptor_channel);
%% Input for pointmatch_task_local_vessel
% This is for generating the descriptor correspondance in z direction 
runlocal = true;
brain = '05_14_58_cube';
inputfolder = sprintf('/data/Vessel/ML_stitching/%s/raw_data',brain);
experimentfolder = sprintf('/data/Vessel/ML_stitching/%s',brain);
descriptorfolder = fullfile(experimentfolder,'stage_2_descriptor_output');
matfolder = fullfile(experimentfolder,'matfiles/');
scopefile = fullfile(matfolder,'scopeloc.mat');
directions = 'Z';
ch = '0';
matchfolder = fullfile(experimentfolder, 'stage_3_point_match_output');
%% Input for pipeline-stitching main_vessel
clc;clear;
runlocal = true;
brain = '05_14_58_cube';
inputfolder = sprintf('/data/Vessel/ML_stitching/%s/raw_data',brain);
experimentfolder = sprintf('/data/Vessel/ML_stitching/%s',brain);
descoutput = fullfile(experimentfolder,'stage_2_descriptor_output');
descriptorfolder = descoutput;
matchoutput = fullfile(experimentfolder, 'stage_3_point_match_output');
matchfolder = matchoutput;
matfolder = fullfile(experimentfolder,'matfiles/');
if ~isfolder(matfolder)
    mkdir(matfolder)
end
scopefile = fullfile(matfolder,'scopeloc.mat');
desc_ch = {'0'};
descriptorfile = fullfile(matfolder,sprintf('descriptors_ch%s.mat',desc_ch{:})); % accumulated descriptor file
runfull = true;
%% Visualization matching 
scatter3(X_(:,1), X_(:,2), X_(:,3))
hold on
scatter3(Y_(:,1) + pixshiftout(1), Y_(:,2) + pixshiftout(2), Y_(:,3) + pixshiftout(3))
%% Analyze descriptor result
[pixel_shift_fft_x, pixel_shift_fft_y, pixel_shift_skel_x, pixel_shift_skel_y, ...
    pixel_shift_edge_x, pixel_shift_edge_y, pixel_shift_stage_y, pixel_shift_stage_x, ...
    avg_dev_stage_x, avg_dev_stage_y, avg_dev_fft_x, avg_dev_fft_y] = deal(zeros(3, numel(paireddescriptor{1})));
[match_rate_fft_x, match_rate_fft_y] = deal(zeros(1, numel(paireddescriptor{1})));
for iter_pair = 1 : numel(paireddescriptor{1})
    tmp_str = paireddescriptor{1}{iter_pair};
    if ~isempty(tmp_str.onx.pixshift_mask_fft)
        pixel_shift_fft_x(:,iter_pair) = tmp_str.onx.pixshift_mask_fft;
        match_rate_fft_x(iter_pair) = tmp_str.onx.match_rate_mask_fft;
    end
    if ~isempty(tmp_str.ony.pixshift_mask_fft)
        pixel_shift_fft_y(:,iter_pair) = tmp_str.ony.pixshift_mask_fft;
        match_rate_fft_y(iter_pair) = tmp_str.ony.match_rate_mask_fft;
    end
    if ~isempty(tmp_str.onx.pixshift_edge)
        pixel_shift_edge_x(:,iter_pair) = tmp_str.onx.pixshift_edge;
    end
    if ~isempty(tmp_str.ony.pixshift_edge)
        pixel_shift_edge_y(:,iter_pair) = tmp_str.ony.pixshift_edge;
    end
    if ~isempty(tmp_str.onx.pixshift_skl)
        pixel_shift_skel_x(:,iter_pair) = tmp_str.onx.pixshift_skl;
    end
    if ~isempty(tmp_str.ony.pixshift_skl)
        pixel_shift_skel_y(:,iter_pair) = tmp_str.ony.pixshift_skl;
    end    
    if isfield(tmp_str.ony, 'pixshift_stage')
        if ~isempty(tmp_str.ony.pixshift_stage)
            pixel_shift_stage_y(:,iter_pair) = tmp_str.ony.pixshift_stage;
        end
    end 
    if isfield(tmp_str.onx, 'pixshift_stage')
        if ~isempty(tmp_str.onx.pixshift_stage)
            pixel_shift_stage_x(:,iter_pair) = tmp_str.onx.pixshift_stage;
        end
    end 
    
    if isfield(tmp_str.onx, 'X') && ~isempty(tmp_str.onx.X)
        tmp_dev = tmp_str.onx.Y - tmp_str.onx.X;
        tmp_dev_stage = tmp_dev - tmp_str.onx.pixshift_stage;
        tmp_dev_fft = tmp_dev - tmp_str.onx.pixshift_mask_fft;
        
        avg_dev_stage_x(:, iter_pair) = mean(tmp_dev_stage, 1);
        avg_dev_fft_x(:, iter_pair) = mean(tmp_dev_fft, 1);
    end
    
    if isfield(tmp_str.ony, 'Y') && ~isempty(tmp_str.ony.X)
        tmp_dev = tmp_str.ony.Y - tmp_str.ony.X;
        tmp_dev_stage = tmp_dev - tmp_str.ony.pixshift_stage;
        tmp_dev_fft = tmp_dev - tmp_str.ony.pixshift_mask_fft;
        avg_dev_stage_y(:, iter_pair) = mean(tmp_dev_stage, 1);
        avg_dev_fft_y(:, iter_pair) = mean(tmp_dev_fft, 1);
    end
end
raw_voxel_size = [0.2969, 0.3758, 1.0000];
avg_dev_stage_x_norm = sqrt(sum(bsxfun(@times, avg_dev_stage_x, raw_voxel_size').^2, 1));
avg_dev_stage_y_norm = sqrt(sum(bsxfun(@times, avg_dev_stage_y, raw_voxel_size').^2, 1));
avg_dev_fft_x_norm = sqrt(sum(bsxfun(@times, avg_dev_fft_x, raw_voxel_size').^2, 1));
avg_dev_fft_y_norm = sqrt(sum(bsxfun(@times, avg_dev_fft_y, raw_voxel_size').^2, 1));
% Deviation from the stage pixel shift vs the fft pixel shift
tmp_str = paireddescriptor{1}{19};
figure;
subplot(2,2,1)
scatter(avg_dev_stage_x_norm, avg_dev_fft_x_norm);
xlabel('Average deviation from stage estimation/\mum');
ylabel('Average deviation from masked fft estimation/\mum');
grid on 
title('Match on X direction')
subplot(2,2,2)
scatter(match_rate_fft_x, avg_dev_fft_x_norm./avg_dev_stage_x_norm)
grid on
xlabel('Cross-correlation');
ylabel('Ratio between masked fft estimation and stage estimation');
title('Match on X direction')
subplot(2,2,3)
scatter(avg_dev_stage_y_norm, avg_dev_fft_y_norm);
xlabel('Average deviation from stage estimation/\mum');
ylabel('Average deviation from masked fft estimation/\mum');
title('Match on Y direction');
grid on 
subplot(2,2,4)
scatter(match_rate_fft_y, avg_dev_fft_y_norm./avg_dev_stage_y_norm)
grid on
xlabel('Cross-correlation');
ylabel('Ratio between masked fft estimation and stage estimation');
title('Match on Y direction');
%% Check descriptor matches
load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
load(fullfile(matfolder,'scopeparams_pertile'))
load(fullfile(matfolder,'regpts'),'regpts')




