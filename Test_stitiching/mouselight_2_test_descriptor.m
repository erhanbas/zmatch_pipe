data_info = load('mouselight_2_ch_0_raw_data_info.mat');
grid_min = min(data_info.stage_grid_xyz);
%%
test_grid_pos = [188, 128, 33];
test_grid_pos_r = test_grid_pos - grid_min + 1;
image_fp = data_info.image_filepath_array{test_grid_pos_r(2), test_grid_pos_r(1), test_grid_pos_r(3)};
im = deployedtiffread(image_fp);
% Light intensity in x, y, z
mean_1 = mean(im, 1);
im_z = squeeze(mean(mean_1, 2));
im_x = squeeze(mean(mean_1, 3));
im_y = squeeze(mean(mean(im, 2), 3));
figure;
subplot(1,3,1);
plot(im_x);
subplot(1,3,2);
plot(im_y);
subplot(1,3,3)
plot(im_z);
implay(im);
%% Intensity Histogram 
histogram(im, 'Normalization', 'cdf')
set(gca,'YScale', 'log');