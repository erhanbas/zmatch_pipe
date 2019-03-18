clc;clear;
% This script is for running the pipeline on specific layer / tiles on the
% Janelia cluster. 
% Taks: 
% 1. Rerun feature extraction, z-matching and x/y matching for section 1541
% - 1544. In those sections, there was a bubble at the middle of the field
% of view and therefore the objects are much dimmer than the normal images.
% The intensity threshold used for thresholding the vessels does not work
% very well with these sections. 
brain  = '2018-08-15';
raw_data_folder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/acquisition/%s', brain);
pipelineoutputfolder = '/nrs/mouselight/pipeline_output/2018-08-15_pipeline_test';
line_fix_folder = fullfile(pipelineoutputfolder, 'stage_1_line_fix_output');
descriptor_folder = fullfile(pipelineoutputfolder, 'stage_2_descriptor_output');
point_match_folder = fullfile(pipelineoutputfolder, 'stage_3_point_match_output');
raw_data_grid = load('mouselight_1_raw_data_info.mat');

experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s_xj_wholebrain', brain);
matfolder = fullfile(experimentfolder,'matfiles/');
directionMap = containers.Map({'-X','-Y','X','Y','-Z','Z'},[ 2, 3, 4, 5, 6, 7]);
scopefile = fullfile(matfolder,'scopeloc.mat');
scopeloc = load(scopefile,'scopeloc');
scopeloc = scopeloc.scopeloc;
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
    
% raw_data_info = load(fullfile(DataManager.SCRIPT_PATH, 'Metadata',  'mouselight_1_raw_data_info.mat'));
configfile = '/groups/mousebrainmicro/home/jix/Documents/GitHub/pipeline-descriptor/configfiles/2018-08-15.cfg';
%% Generate descriptor
compute_tile_idx = find(raw_data_grid.stage_grid_xyz(:,3) >= 1540 & ...
    raw_data_grid.stage_grid_xyz(:,3) <= 1544);
num_file = numel(compute_tile_idx);

num_core = feature('numcores');
poolobj = gcp('nocreate');
delete(poolobj);
parpool(num_core);
parfor file_idx = 1 : num_file
    tmp_tile_idx = compute_tile_idx(file_idx);
    raw_image_fp = raw_data_grid.image_filepath{tmp_tile_idx};
    % Convert to line fix data
    image_fp = strrep(raw_image_fp, raw_data_folder, line_fix_folder);
    output_fp = strrep(raw_image_fp, raw_data_folder, descriptor_folder);    
    output_fp = strrep(output_fp, 'ngc.0.tif', 'desc.0.mat');
    output_folder = fileparts(output_fp);
    fprintf('Processing %s (%d/%d) \n', image_fp, file_idx, num_file);
    if ~isfolder(output_folder)
        mkdir(output_folder);
    end
    try
        exit_code = vesselDescriptor(image_fp, output_fp, configfile);
    catch ME
        sprintf('Fail to process tiel %s\n', image_fp);
        sprintf('Error message: %s\n', ME.identifier);
    end
end
disp('Finish computing descriptor');
%% Descriptor match in z-direction 
match_z_idx = find(scopeloc.gridix(:,3) >= 1540 & ...
    scopeloc.gridix(:,3) <= 1544);
num_z_match_tile = numel(match_z_idx);
num_core = feature('numcores');
poolobj = gcp('nocreate');
delete(poolobj);
parpool(round(num_core * 1));
directions = 'Z';
ch = {'0'};
maxnumofdesc = 1e4;
parfor iter_tile = 1 : num_z_match_tile
    fprintf('Processing tile %d/%d\n', iter_tile, num_z_match_tile);
    tmp_tile_idx = match_z_idx(iter_tile);
%     tmp_tile_idx = test_grid_idx;
    tile1 = fullfile(descriptor_folder, scopeloc.relativepaths{tmp_tile_idx});
    acqusitionfolder1 = fileparts(scopeloc.filepath{tmp_tile_idx});
    iineig = neighbors(tmp_tile_idx, directionMap(directions));
    outfold =fullfile(point_match_folder,scopeloc.relativepaths{tmp_tile_idx});
    if ~isnan(iineig)
        tile2 = fullfile(descriptor_folder,scopeloc.relativepaths{iineig});
        acqusitionfolder2 = fileparts(scopeloc.filepath{iineig});
    else
%         fprintf('Tile %d does not has a neighbor below\n', tmp_tile_idx);
%         parfor_progress;
        continue;
    end
%     if isfile(fullfile(outfold, 'match-Z.mat'))
% %         fprintf('Tile %d has been processed\n', tmp_tile_idx);
% %         parfor_progress;
%         continue;        
%     end
    try
        pointmatch_vessel(tile1,tile2,acqusitionfolder1,acqusitionfolder2,outfold,[0,0,0],ch,maxnumofdesc,0);
    catch ME
        sprintf('Fail to process tiel %s\n', tile1);
        sprintf('Error message: %s\n', ME.identifier);
    end
%     parfor_progress;
end
disp('Finished Z-matching');