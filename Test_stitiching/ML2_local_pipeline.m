data_info = load('mouselight_2_ch_0_line_fix_data_info.mat');
raw_data_folder = '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2019-01-24';
dataset_folder = '/nrs/mouselight/pipeline_output/2019-01-24';
line_fix_folder_name = 'stage_1_line_fix_output';
line_fix_folder = fullfile('/nrs/mouselight/pipeline_output/2019-01-24', line_fix_folder_name);
descriptor_folder_name = 'stage_2_descriptor_output';
descriptor_folder = fullfile(dataset_folder, descriptor_folder_name);

local_folder = '/data/Vessel/ML_stitching/2019-01-24';
local_line_fix_folder = fullfile(local_folder, line_fix_folder_name);
local_descriptor_folder = fullfile(local_folder, [descriptor_folder_name, '_v2']);

config_file = '/home/dklab/Documents/Github/MouseLight_repos/pipeline-descriptor/configfiles/2018-08-15.cfg';
%% Specify the test region 
test_region_grid_coordinate = [194, 129, 26];
test_region_list_ind = find(~any(data_info.stage_grid_xyz - test_region_grid_coordinate, 2));
test_region_grid_pos_sub = data_info.grid_pos_sub(test_region_list_ind, :);
num_neighbor = 1;
% Here should use meshgrid instead of ndgrid
[region_pos_sub_1, region_pos_sub_2, region_pos_sub_3] = meshgrid((test_region_grid_pos_sub(1)-num_neighbor):(test_region_grid_pos_sub(1)+num_neighbor), ...
    (test_region_grid_pos_sub(2)-num_neighbor):(test_region_grid_pos_sub(2)+num_neighbor), ...
    (test_region_grid_pos_sub(3)-num_neighbor):(test_region_grid_pos_sub(3)+num_neighbor));
region_pos_sub_1 = region_pos_sub_1(:);
region_pos_sub_2 = region_pos_sub_2(:);
region_pos_sub_3 = region_pos_sub_3(:);
region_pos_ind = sub2ind(data_info.grid_size, region_pos_sub_1, region_pos_sub_2, region_pos_sub_3);
region_linear_idx = data_info.linear_idx_array(region_pos_ind);

test_region_file_list = data_info.image_filepath_array(region_pos_ind);
region_w_file_Q = ~cellfun(@isempty, test_region_file_list);
num_block = numel(test_region_file_list);
num_file = nnz(region_w_file_Q);
%% Copy files to local machine
for iter_file = 1 : num_block
    tmp_fp = test_region_file_list{iter_file};
    if isempty(tmp_fp)
        continue;
    end
    tmp_folder = fileparts(tmp_fp);
    tmp_folder = strrep(tmp_folder, raw_data_folder, line_fix_folder);
    if isfolder(tmp_folder)
       tmp_target_subfolder = strrep(tmp_folder, line_fix_folder, local_line_fix_folder);
       if ~isfolder(tmp_target_subfolder)
           mkdir(tmp_target_subfolder);
       end
       % Copy the entire folder
       fun_copy_by_rsync(tmp_folder, tmp_target_subfolder, true, false);        
    end    
end
%%




%% Generate descriptors
parfor (file_idx = 1 : num_block, 9)
    tmp_fp = test_region_file_list{file_idx};
    if isempty(tmp_fp)
        continue;
    end
    image_fp = strrep(tmp_fp, raw_data_folder, local_line_fix_folder);
    
    [output_folder, output_file_name, ~] = fileparts(image_fp);
    output_folder = strrep(output_folder, local_line_fix_folder, local_descriptor_folder);
    output_file_name = sprintf('%s.mat', strrep(output_file_name, 'ngc', 'desc'));
    output_fp = fullfile(output_folder, output_file_name);
    
    fprintf('Processing %s (%d/%d) \n', image_fp, file_idx, num_block);
    if ~isfolder(output_folder)
        mkdir(output_folder);
    end
    tmp_tic = tic;
    exit_code = vesselDescriptor_v2(image_fp, output_fp, config_file);
%     exit_code = skelDescriptor(image_fp, output_fp, configfile);
    fprintf('Finish processing tile %d / %d. Elapsed time was %f seconds.\n', file_idx, ...
        num_block, toc(tmp_tic));
end
disp('Finish computing descriptor');
%% Generate z-matching 
% Generate the scope location file - or use my data_info structure? 
scope_acquisition_info_fp = fullfile(local_line_fix_folder, 'scopeacquisitionlist.txt');
if ~isfile(scope_acquisition_info_fp)
    scope_files = dir(fullfile(local_line_fix_folder, ['**', filesep, '*.acquisition']));
    scope_file_names = cellfun(@fullfile, {scope_files.folder}, {scope_files.name}, 'UniformOutput', false);
    % Output scope acquisition list ( required by pointmatch_task
    fid = fopen(scope_acquisition_info_fp, 'w');
    for idx = 1 : numel(scope_file_names)
        fprintf(fid, '%s\n', scope_file_names{idx});
    end
    fclose(fid);
end
[scopeloc] = getScopeCoordinates(local_line_fix_folder, 1);

%% Create scopeloc by data_info
new_scopeloc = struct;
new_scopeloc.gridix = nan(num_file, 4);
new_scopeloc.gridix(:, 1:3) = data_info.stage_grid_xyz(region_linear_idx, :);
new_scopeloc.loc = data_info.stage_yzx_um([2,1,3],region_pos_ind)'; 

new_scopeloc.filepath = test_region_file_list(region_w_file_Q);
new_scopeloc.relativepaths = cellfun(@(x) fileparts(erase(x, raw_data_folder)), new_scopeloc.filepath, 'UniformOutput', false);








%% Generate x/y matching


%% Visualization of descriptor matching 
%% Debuging
% DataManager = FileManager;
% vis_mask = uint8(Io_mask);
% vis_mask(kept_skl_ind) = 2;
% DataManager.visualize_itksnap(Io, vis_mask);