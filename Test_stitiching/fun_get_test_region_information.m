function info = fun_get_test_region_information(data_info, region_grid_ind, raw_data_root, test_data_root)


line_fix_folder_name = 'stage_1_line_fix_output';
descriptor_folder_name = 'stage_2_descriptor_output';


region_linear_idx = data_info.linear_idx_array(region_grid_ind);

info = struct;
info.raw_data_file_list = data_info.image_filepath_array(region_grid_ind);
info.raw_data_file_exist_Q = ~cellfun(@isempty, info.raw_data_file_list);
info.num_file = nnz(info.raw_data_file_exist_Q);
for iter_file = 1 : info.num_block





end