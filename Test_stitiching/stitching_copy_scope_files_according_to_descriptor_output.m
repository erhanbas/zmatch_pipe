raw_data_root_folder = '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2018-08-15';
pipeline_output_folder = '/nrs/mouselight/pipeline_output/2018-08-15_pipeline_test';
pipeline_scope_folder = fullfile(pipeline_output_folder, 'raw_data_info');
descriptor_output_folder = fullfile(pipeline_output_folder, 'stage_2_descriptor_output'); 
descriptor_file_str = dir(fullfile(descriptor_output_folder, '**', '*.mat'));
for iter_file = 1 : numel(descriptor_file_str)
    fprintf('Processing file %d\n', iter_file);
    tmp_descriptor_fp = fullfile(descriptor_file_str(iter_file).folder, ...
        descriptor_file_str(iter_file).name);
    tmp_descriptor_folder = fileparts(tmp_descriptor_fp);
    tmp_dest_tile_folder = strrep(tmp_descriptor_fp, 'stage_2_descriptor_output', 'raw_data_info');
    tmp_dest_tile_folder = fileparts(tmp_dest_tile_folder);
    if ~isfolder(tmp_dest_tile_folder)
        mkdir(tmp_dest_tile_folder)
    end
    tmp_source_fp = strrep(tmp_descriptor_folder, descriptor_output_folder, raw_data_root_folder);
    tmp_acquisition_fp = dir(fullfile(tmp_source_fp, '*acquisition*'));
    tmp_dest_acquisition_fp = fullfile(tmp_dest_tile_folder, tmp_acquisition_fp.name);
    tmp_source_acquisition_fp = fullfile(tmp_acquisition_fp.folder, tmp_acquisition_fp.name);
    
    tmp_scope_fp_str = dir(fullfile(tmp_source_fp, '*.microscope*'));
    tmp_dest_scope_fp = fullfile(tmp_dest_tile_folder, tmp_scope_fp_str.name);
    tmp_source_scope_fp = fullfile(tmp_scope_fp_str.folder, tmp_scope_fp_str.name);
    % Copy file
    
    copyfile(tmp_source_acquisition_fp, tmp_dest_acquisition_fp);
    copyfile(tmp_source_scope_fp, tmp_dest_scope_fp);
end