% clc;clear;
data_info = load('mouselight_2_ch_0_line_fix_data_info.mat');

raw_data_folder = '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2019-01-24';
dataset_folder = '/nrs/mouselight/pipeline_output/2019-01-24';
line_fix_folder = '/nrs/mouselight/pipeline_output/2019-01-24/stage_1_line_fix_output';
descriptor_folder = fullfile(dataset_folder, 'stage_2_descriptor_output');
num_file = data_info.num_tile;
% dataset_folder = fullfile(Vessel_data_root_path, 'pipeline_test', '06_19_58_cube');
% raw_data_folder = fullfile(dataset_folder, 'raw_data');
% input_file_path = dir(fullfile(raw_data_folder, ['**', filesep, '*.tif']));
% scope_files = dir(fullfile(raw_data_folder, ['**', filesep, '*.acquisition']));
% scope_file_names = cellfun(@fullfile, {scope_files.folder}, {scope_files.name}, 'UniformOutput', false);
% % Output scope acquisition list ( required by pointmatch_task
% fid = fopen(fullfile(raw_data_folder, 'scopeacquisitionlist.txt'), 'w');
% for idx = 1 : numel(scope_file_names)
%     fprintf(fid, '%s\n', scope_file_names{idx});
% end
% fclose(fid);
% num_file = numel(input_file_path);
configfile = '/groups/mousebrainmicro/home/jix/Documents/GitHub/pipeline-descriptor/configfiles/2018-08-15.cfg';
%% Debug stitching
debug_grid_xyz = [193 129 27];
debug_list_ind = find(~any(data_info.stage_grid_xyz - debug_grid_xyz, 2));

image_fp = strrep(data_info.image_filepath{debug_list_ind}, raw_data_folder, line_fix_folder);
[output_folder, output_file_name, ~] = fileparts(image_fp);
output_folder = strrep(output_folder, line_fix_folder, descriptor_folder);
output_file_name = sprintf('%s.mat', strrep(output_file_name, 'ngc', 'desc'));
output_fp = fullfile(output_folder, output_file_name);

tmp_desc = load(output_fp);
% Matched result
debug_match_Z_ds = regpts{debug_list_ind};
debug_match_Z = featmap(debug_list_ind).Z;
% Image
tmp_im = deployedtiffread(image_fp);
%% Generate descriptor
% image_fp = fullfile('/nrs/mouselight/pipeline_output/2019-01-24/stage_1_line_fix_output',...
%     '2019-02-18', '01', '01989', '01989-ngc.0.tif');
tmp_im = deployedtiffread(image_fp);
file_idx = debug_list_ind;
for file_idx = 1 : num_file
    image_fp = strrep(data_info.image_filepath{file_idx}, raw_data_folder, line_fix_folder);
    [output_folder, output_file_name, ~] = fileparts(image_fp);
    output_folder = strrep(output_folder, line_fix_folder, descriptor_folder);
    output_file_name = sprintf('%s.mat', strrep(output_file_name, 'ngc', 'desc'));
    output_fp = fullfile(output_folder, output_file_name);
    
    fprintf('Processing %s (%d/%d) \n', image_fp, file_idx, num_file);
    if isfile(output_fp)
        tmp_str = load(output_fp, 'record');
        if ~isempty(tmp_str.record)
            continue;
        end
    end
    if ~isfolder(output_folder)
        mkdir(output_folder);
    end
    tic
    exit_code = vesselDescriptor(image_fp, output_fp, configfile);
%     exit_code = skelDescriptor(image_fp, output_fp, configfile);
    toc
end
disp('Finish computing descriptor');
%% Point Match - between two tiles
matching_folder = fullfile(dataset_folder, 'stage_3_point_match_output');
tile_1_grid_xyz = [193, 129, 27];
tile_2_grid_xyz = [193, 129, 28];
tile_idx_1 = find(~any(data_info.stage_grid_xyz - tile_1_grid_xyz, 2));
tile_idx_2 = find(~any(data_info.stage_grid_xyz - tile_2_grid_xyz, 2));

tile_descriptor_fp_1 = strrep(strrep(data_info.image_filepath{tile_idx_1}, raw_data_folder, descriptor_folder), 'ngc.0.tif', 'desc.0.mat');
tile_descriptor_fp_2 = strrep(strrep(data_info.image_filepath{tile_idx_2}, raw_data_folder, descriptor_folder), 'ngc.0.tif', 'desc.0.mat');
tile_acqusition_fd_1 = fileparts(data_info.image_filepath{tile_idx_1});
tile_acqusition_fd_2 = fileparts(data_info.image_filepath{tile_idx_2});
output_folder = strrep(fileparts(tile_descriptor_fp_1), descriptor_folder, matching_folder);
if ~isfolder(output_folder)
    mkdir(output_folder);
end
descriptor_channel = {'0'};
pixshift = [0,0,0];
[~] = pointmatch_vessel(tile_descriptor_fp_1, tile_descriptor_fp_2, tile_acqusition_fd_1, tile_acqusition_fd_2, ...
    output_folder,pixshift, descriptor_channel, 1e4, 0);
%% Input for pointmatch_task_local_vessel
clc;clear
runlocal = true;
brain = '2018-08-15';
% inputfolder is for getting the scope location 
inputfolder = '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2018-08-15';
% Experiment folder is for writing the control points, matfiles
experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s_xj_wholebrain', brain);
% Pipeline output folder
pipelineoutputfolder = '/nrs/mouselight/pipeline_output/2018-08-15_pipeline_test';
descriptorfolder = fullfile(pipelineoutputfolder, 'stage_2_descriptor_output');
matchfolder = fullfile(pipelineoutputfolder, 'stage_3_point_match_output');

matfolder = fullfile(experimentfolder,'matfiles/');
scopefile = fullfile(matfolder,'scopeloc.mat');
directions = 'Z';
ch = '0';

%% Input for pipeline-stitching main.m
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






