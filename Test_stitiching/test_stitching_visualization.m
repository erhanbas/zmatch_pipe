DataManager = FileManager;
vis_mask = uint8(Io_mask);
vis_mask(skel) = 3;
vis_mask(kept_skl_ind) = 2;

DataManager.visualize_itksnap(Io, vis_mask)
%% Matched skeleton voxels
figure;
scatter3(X_skel(:,1), X_skel(:,2), X_skel(:,3))
hold on 
scatter3(Y_skel(:,1), Y_skel(:,2), Y_skel(:,3))
legend('Tile 1', 'Tile 2', 'Tile 1 matched', 'Tile 2 matched');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(2)
%% Input surface voxels
%% Matched surface voxels
figure;
scatter3(X_edge(:,1), X_edge(:,2), X_edge(:,3))
hold on 
scatter3(Y_edge(:,1), Y_edge(:,2), Y_edge(:,3))
hold on
scatter3(X_skel(:,1), X_skel(:,2), X_skel(:,3))
hold on 
scatter3(Y_skel(:,1), Y_skel(:,2), Y_skel(:,3))
legend('Tile 1 edge', 'Tile 2 edge', 'Tile 1 skeleton', 'Tile 2 skeleton');
xlabel('X');
ylabel('Y');
zlabel('Z');

%% Plot selected points
figure;
subplot(1,3,1)
scatter3(des_fixed(:,1), des_fixed(:,2), des_fixed(:,3))
hold on 
scatter3(des_moving_ori(:,1), des_moving_ori(:,2), des_moving_ori(:,3))
legend('Tile 1', 'Tile 2', 'Tile 1 matched', 'Tile 2 matched');
xlabel('X');
ylabel('Y');
zlabel('Z');
% view(2)
subplot(1,3,2)
scatter3(des_1_sub(:,1), des_1_sub(:,2), des_1_sub(:,3))
hold on 
scatter3(des_2_sub_shift(:,1), des_2_sub_shift(:,2), des_2_sub_shift(:,3))
legend('Tile 1', 'Tile 2', 'Tile 1 matched', 'Tile 2 matched');
xlabel('X');
ylabel('Y');
zlabel('Z');
% view(2)
subplot(1,3,3)
scatter3(X_(:,1), X_(:,2), X_(:,3))
hold on 
scatter3(Y_(:,1), Y_(:,2), Y_(:,3))
legend('Tile 1 matched', 'Tile 2 matched');
xlabel('X');
ylabel('Y');
zlabel('Z');
% view(2)
%% Plot all input point ( shifted ) 
figure;
scatter3(des_1_sub(:,1), des_1_sub(:,2), des_1_sub(:,3))
hold on 
scatter3(des_2_sub_shift(:,1), des_2_sub_shift(:,2), des_2_sub_shift(:,3))
%% Edge downsampled shifted
figure;
scatter3(desc_1_sub_ds(:,1), desc_1_sub_ds(:,2), desc_1_sub_ds(:,3))
hold on 
scatter3(desc_2_sub_ds(:,1), desc_2_sub_ds(:,2), desc_2_sub_ds(:,3))
%%
figure;
scatter3(des_center(:,1), des_center(:,2), des_center(:,3))
hold on 
scatter3(des_cadj_ori(:,1), des_cadj_ori(:,2), des_cadj_ori(:,3))
legend('Tile 1', 'Tile 2');
xlabel('X');
ylabel('Y');
zlabel('Z');
%% Plot matched points
figure;
% subplot(1,2,1)
scatter3(X_(:,1), X_(:,2), X_(:,3))
hold on 
% scatter3(Y_(:,1), Y_(:,2), Y_(:,3))
% hold on 
scatter3(Y_(:,1) + pixshift(1), Y_(:,2) + pixshift(2), Y_(:,3) + pixshift(3))
% title('Matched edge voxels')
% legend('Tile 1', 'Tile 2', 'Tile 2 moved');
xlabel('X');
ylabel('Y');
zlabel('Z');
%%
figure;
scatter3(des_center(:,1), des_center(:,2), des_center(:,3))
hold on 
scatter3(des_cadj_ori(:,1), des_cadj_ori(:,2), des_cadj_ori(:,3))
legend('Tile 1', 'Tile 2');
xlabel('X');
ylabel('Y');
zlabel('Z');
% view(2)
% scatter3(tY_(:,1), tY_(:,2), tY_(:,3))
%%
% figure;
% scatter3(desc1(:,1), desc1(:,2), desc1(:,3))
% hold on 
% scatter3(desc2(:,1), desc2(:,2), desc2(:,3))
% legend('Tile 1', 'Tile 2');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
%%
figure;
scatter3(desc1_skel(:,1), desc1_skel(:,2), desc1_skel(:,3))
hold on 
scatter3(desc2_skel(:,1), desc2_skel(:,2), desc2_skel(:,3))
legend('Tile 1', 'Tile 2');
xlabel('X');
ylabel('Y');
zlabel('Z');
% view(2)
%%
figure;
scatter3(desc1(:,1), desc1(:,2), desc1(:,3))
hold on 
scatter3(desc2(:,1), desc2(:,2), desc2(:,3))
legend('Tile 1', 'Tile 2');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(2)
%%
figure;
scatter3(X_sampled(:,1), X_sampled(:,2), X_sampled(:,3))
hold on 
scatter3(Y_sampled(:,1), Y_sampled(:,2), Y_sampled(:,3))
legend('Tile 1', 'Tile 2');
xlabel('X');
ylabel('Y');
zlabel('Z');
% view(2)
%%
figure;
scatter3(descriptor_1_sub(1:10:end,1), descriptor_1_sub(1:10:end,2), descriptor_1_sub(1:10:end,3))
hold on 
scatter3(desc_2_sub_shifted(1:10:end,1), desc_2_sub_shifted(1:10:end,2), desc_2_sub_shifted(1:10:end,3))
legend('Tile 1', 'Tile 2');
xlabel('X');
ylabel('Y');
zlabel('Z');
%%

figure;
scatter3(X_stable(1:10:end,1), X_stable(1:10:end,2), X_stable(1:10:end,3))
hold on 
scatter3(Y_stable(1:10:end,1), Y_stable(1:10:end,2), Y_stable(1:10:end,3))
legend('Tile 1', 'Tile 2');
xlabel('X');
ylabel('Y');
zlabel('Z');
%%
figure;
scatter3(allX(1:10:end,1), allX(1:10:end,2), allX(1:10:end,3))
hold on 
scatter3(allY(1:10:end,1), allY(1:10:end,2), allY(1:10:end,3))
legend('Tile 1', 'Tile 2');
xlabel('X');
ylabel('Y');
zlabel('Z');
%%
tmpX_edge = tmp_str.onx.X_edge;
tmpY_edge = tmp_str.onx.Y_edge;
tmpX_skel = tmp_str.onx.X_skl;
tmpY_skel = tmp_str.onx.Y_skl;

figure;
scatter3(tmpX_edge(:,1), tmpX_edge(:,2), tmpX_edge(:,3));
hold on 
scatter3(tmpY_edge(:,1), tmpY_edge(:,2), tmpY_edge(:,3));
% figure
hold on
scatter3(tmpX_skel(:,1), tmpX_skel(:,2), tmpX_skel(:,3));
hold on 
scatter3(tmpY_skel(:,1), tmpY_skel(:,2), tmpY_skel(:,3));
legend('X edge', 'Y edge', 'X skeleton', 'Y skeleton')
%% Visualize the strong gradient edge voxels
for iter_vis = 202 : numel(paireddescriptor{1})
    pause(3);
    %%
    iter_vis = 282;
    vis_tile_str = paireddescriptor{1}{iter_vis };
    vis_tile_str_x = vis_tile_str.onx;
    vis_tile_str_y = vis_tile_str.ony;
%     fprintf('================================\n');
    fprintf('********Pair %d ****************\n', iter_vis);
    disp('Match in X direction');
    if ~isempty(vis_tile_str_x.pixshift_stage)
        fprintf('Pixel shift stage: (%d, %d, %d)\n', vis_tile_str_x.pixshift_stage);
    end
    if ~isempty(vis_tile_str_x.pixshift_mask_fft)
        fprintf('Pixel shift Masked_fft: (%d, %d, %d)\n', vis_tile_str_x.pixshift_mask_fft);
        fprintf('Score Masked_fft: %f\n', vis_tile_str_x.match_rate_mask_fft);
    end
    if ~isempty(vis_tile_str_x.pixshift_skl)
        fprintf('Pixel shift skel: (%d, %d, %d)\n', vis_tile_str_x.pixshift_skl);
        fprintf('Score skel: %f\n', vis_tile_str_x.match_rate_skl);
    end
    if ~isempty(vis_tile_str_x.pixshift_edge)
        fprintf('Pixel shift edge: (%d, %d, %d)\n', vis_tile_str_x.pixshift_edge);
        fprintf('Score edge: %f\n', vis_tile_str_x.match_rate_edge);
    end
    fprintf('================================\n');
    disp('Match in Y direction');
    if ~isempty(vis_tile_str_y.pixshift_stage)
        fprintf('Pixel shift stage: (%d, %d, %d)\n', vis_tile_str_y.pixshift_stage);
    end
    if ~isempty(vis_tile_str_y.pixshift_mask_fft)
        fprintf('Pixel shift Masked_fft: (%d, %d, %d)\n', vis_tile_str_y.pixshift_mask_fft);
        fprintf('Score Masked_fft: %f\n', vis_tile_str_y.match_rate_mask_fft);
    end
    if ~isempty(vis_tile_str_y.pixshift_skl)
        fprintf('Pixel shift skel: (%d, %d, %d)\n', vis_tile_str_y.pixshift_skl);
        fprintf('Score skel: %f\n', vis_tile_str_y.match_rate_skl);
    end
    if ~isempty(vis_tile_str_y.pixshift_edge)
        fprintf('Pixel shift edge: (%d, %d, %d)\n', vis_tile_str_y.pixshift_edge);
        fprintf('Score edge: %f\n', vis_tile_str_y.match_rate_edge);
    end
    % figure;
    subplot(2,2,1)
    cla;
    if ~isempty(vis_tile_str_x.X_skl)
        scatter3(vis_tile_str_x.X_skl(:,1), vis_tile_str_x.X_skl(:,2), vis_tile_str_x.X_skl(:,3));
        hold on 
        scatter3(vis_tile_str_x.Y_skl(:,1), vis_tile_str_x.Y_skl(:,2), vis_tile_str_x.Y_skl(:,3));
%         figure
%         scatter3(vis_tile_str_x.X_skl(:,1), vis_tile_str_x.X_skl(:,2), vis_tile_str_x.X_skl(:,3));
%         hold on 
%         scatter3(vis_tile_str_x.Y_skl(:,1) - vis_tile_str_x.X_skl(:,1),...
%             vis_tile_str_x.Y_skl(:,2) - vis_tile_str_x.X_skl(:,2),...
%             vis_tile_str_x.Y_skl(:,3) - vis_tile_str_x.X_skl(:,3));
        xlim([1, 1024]);
        ylim([1, 1536]);
        zlim([1, 251]);
        daspect([1,1,1])
        xlabel('X'); ylabel('Y'); zlabel('Z');
    end
    if ~isempty(vis_tile_str_x.X_edge)
        hold on
        scatter3(vis_tile_str_x.X_edge(:,1), vis_tile_str_x.X_edge(:,2), vis_tile_str_x.X_edge(:,3));
        hold on 
        scatter3(vis_tile_str_x.Y_edge(:,1), vis_tile_str_x.Y_edge(:,2), vis_tile_str_x.Y_edge(:,3));
    %     legend('Tile 1 edge', 'Tile 2 edge');
        xlim([1, 1024]);
        ylim([1, 1536]);
        zlim([1, 251]);
        daspect([1,1,1])
        hold off
    end
    legend('Tile 1 skel', 'Tile 2 skel', 'Tile 1 edge', 'Tile 2 edge');
    subplot(2,2,2)
    cla;
    if ~isempty(vis_tile_str_y.X_skl)
        scatter3(vis_tile_str_y.X_skl(:,1), vis_tile_str_y.X_skl(:,2), vis_tile_str_y.X_skl(:,3));
        hold on 
        scatter3(vis_tile_str_y.Y_skl(:,1), vis_tile_str_y.Y_skl(:,2), vis_tile_str_y.Y_skl(:,3));
    %     legend('Tile 1 skel', 'Tile 2 skel');
        xlim([1, 1024]);
        ylim([1, 1536]);
        zlim([1, 251]);
        daspect([1,1,1]);
        xlabel('X'); ylabel('Y'); zlabel('Z');
    end
    if ~isempty(vis_tile_str_y.X_edge)
        hold on
        scatter3(vis_tile_str_y.X_edge(:,1), vis_tile_str_y.X_edge(:,2), vis_tile_str_y.X_edge(:,3));
        hold on 
        scatter3(vis_tile_str_y.Y_edge(:,1), vis_tile_str_y.Y_edge(:,2), vis_tile_str_y.Y_edge(:,3));
    %     legend('Tile 1 edge', 'Tile 2 edge');
        xlim([1, 1024]);
        ylim([1, 1536]);
        zlim([1, 251]);
        daspect([1,1,1])
        hold off
    end
    legend('Tile 1 skel', 'Tile 2 skel', 'Tile 1 edge', 'Tile 2 edge');
       
    subplot(2,2,3)
    cla;
    if ~isempty(vis_tile_str_x.X_skl)
        scatter3(vis_tile_str_x.Y_skl(:,1) - vis_tile_str_x.X_skl(:,1),...
            vis_tile_str_x.Y_skl(:,2) - vis_tile_str_x.X_skl(:,2),...
            vis_tile_str_x.Y_skl(:,3) - vis_tile_str_x.X_skl(:,3));
        daspect([1,1,1])
        xlabel('X'); ylabel('Y'); zlabel('Z');
    end
    if ~isempty(vis_tile_str_x.X_edge)
        hold on
        scatter3(vis_tile_str_x.Y_edge(:,1) - vis_tile_str_x.X_edge(:,1),...
            vis_tile_str_x.Y_edge(:,2) - vis_tile_str_x.X_edge(:,2),...
            vis_tile_str_x.Y_edge(:,3) - vis_tile_str_x.X_edge(:,3));
        daspect([1,1,1])
        hold off
    end
    subplot(2,2,4)
    cla;
    if ~isempty(vis_tile_str_y.X_skl)
        scatter3(vis_tile_str_y.Y_skl(:,1) - vis_tile_str_y.X_skl(:,1),...
            vis_tile_str_y.Y_skl(:,2) - vis_tile_str_y.X_skl(:,2),...
            vis_tile_str_y.Y_skl(:,3) - vis_tile_str_y.X_skl(:,3));
        daspect([1,1,1])
        xlabel('X'); ylabel('Y'); zlabel('Z');
    end
    if ~isempty(vis_tile_str_y.X_edge)
        hold on
        scatter3(vis_tile_str_y.Y_edge(:,1) - vis_tile_str_y.X_edge(:,1),...
            vis_tile_str_y.Y_edge(:,2) - vis_tile_str_y.X_edge(:,2),...
            vis_tile_str_y.Y_edge(:,3) - vis_tile_str_y.X_edge(:,3));
        daspect([1,1,1])
        hold off
    end
    %%
end